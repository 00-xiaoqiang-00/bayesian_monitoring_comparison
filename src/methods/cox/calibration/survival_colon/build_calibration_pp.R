# Generate the two full PP-path inputs for Cox survival_colon parameter selection.
# Run through workflows/calibration/run_cox_survival_colon.R with external output.
if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Use the registered workflow.")
library(survival)
library(brms)
library(posterior)
library(parallel)
if (requireNamespace("RcppParallel",quietly=TRUE)) RcppParallel::setThreadOptions(numThreads=1)
output <- Sys.getenv("V33_OUTPUT_DIR")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
setwd(output)
source(file.path(Sys.getenv("V33_CODE_ROOT"),"src/methods/cox/calibration/survival_colon/data_simulation_final.R"),local=TRUE)
set.seed(20260719)
trial_settings <- survival_colon_trial_settings()
colon_observed <- as.data.frame(data.simulation(N.sim = 1, N.max = trial_settings$N.max, cens_upper = trial_settings$cens_upper, 
    scenario = "H1_observed")[, , 1])
raw_data <- transform(colon_observed, arm = group + 1, death = death_t, recurrence = recurrence_t, censor = censor_t, delta1 = as.integer(recurrence_t <= 
    pmin(death_t, censor_t)), delta2 = as.integer(death_t <= censor_t))
cat(sprintf("成功读取原始数据！总样本量: %d 人 (Arm 1/对照组: %d 人, Arm 2/试验组: %d 人)\n", nrow(raw_data), 
    sum(raw_data$arm == 1), sum(raw_data$arm == 2)))
death_true <- ifelse(is.na(raw_data$death), raw_data$censor, raw_data$death)
recurrence_true <- ifelse(is.na(raw_data$recurrence), raw_data$censor, raw_data$recurrence)
raw_data$composite_time <- pmin(death_true, recurrence_true, raw_data$censor)
raw_data$composite_status <- as.integer(raw_data$delta1 == 1 | raw_data$delta2 == 1)
raw_data$group <- raw_data$arm - 1
N.max <- trial_settings$N.max
interim_fracs <- trial_settings$interim_fracs
recruit.int <- trial_settings$recruit.int
bayes_iter <- 2000
bayes_warmup <- 1000
cache_tag <- "survival_colon_v1_H0_H1calibration_H1observed_cox_iter2000_warmup1000"
get_valid_group_draws <- function(fit) {
    if (inherits(fit, "try-error")) 
        return(NULL)
    draws <- try(posterior::as_draws_array(fit, variable = "b_group"), silent = TRUE)
    if (inherits(draws, "try-error")) 
        return(NULL)
    diagnostics <- posterior::summarise_draws(draws)
    if (any(!is.finite(diagnostics$rhat)) || any(diagnostics$rhat > 1.01) || any(diagnostics$ess_bulk < 400) || any(diagnostics$ess_tail < 
        400)) {
        return(NULL)
    }
    posterior::as_draws_matrix(draws)[, "b_group"]
}
set.seed(42)
raw_data <- raw_data[sample(nrow(raw_data)), ]
raw_data <- raw_data[seq_len(N.max), ]
entry_times <- recruit.int * (0:(N.max - 1))
cat("\n【系统优化】正在预编译基础贝叶斯模型 (以极大提速后续分析)...\n")
df_dummy <- data.frame(time = runif(100, 0.5, 40), status = sample(c(0, 1), 100, replace = TRUE), group = sample(c(0, 1), 
    100, replace = TRUE))
base_model <- brm(formula = time | cens(1 - status) ~ group, data = df_dummy, family = brmsfamily("cox"), prior = prior(normal(0, 
    1.5), class = "b", coef = "group"), chains = 1, cores = 1, iter = 10, warmup = 5, seed = 12345, refresh = 0, silent = 2)
cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if(nzchar(cores_env)) as.integer(cores_env) else 30L
N.sim_realcase <- 10000L
generate_calibration_scenario <- function(label, scenario) {
    scenario_tag <- scenario
    cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
    use_cache <- FALSE
    if (file.exists(cache_file)) {
        stats_list_temp <- readRDS(cache_file)
        if (length(stats_list_temp) == N.sim_realcase) {
            cat(sprintf("\n【系统提示】检测到匹配当前模拟数的 [%s] 缓存，直接读取以节省计算时间...\n", 
                label))
            stats_list <- stats_list_temp
            use_cache <- TRUE
        }
    }
    if (!use_cache) {
        cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
        set.seed(777)
        sim_array <- data.simulation(N.sim = N.sim_realcase, N.max = N.max, cens_upper = trial_settings$cens_upper, scenario = scenario)
        stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
            df_j <- as.data.frame(sim_array[, , j])
            time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
            final_censor_t <- pmin(df_j$censor_t, survival_colon_admin_followup(entry_times, trial_settings$final_analysis_time))
            df_j$time <- pmin(time_to_first, final_censor_t)
            df_j$status <- as.integer(time_to_first <= final_censor_t)
            interim_PP <- numeric(length(interim_fracs))
            original_censor_t <- df_j$censor_t
            simulation_ok <- TRUE
            for (k in seq_along(interim_fracs)) {
                n_k <- floor(N.max * interim_fracs[k])
                Time.current <- recruit.int * (n_k - 1)
                obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
                censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
                time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
                status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
                valid_rows <- censor_t_interim > 0
                df_interim <- data.frame(time = time_interim[valid_rows], status = status_interim[valid_rows], group = df_j$group[1:n_k][valid_rows])
                fit_bayes_interim <- try(update(base_model, newdata = df_interim, chains = 4, cores = 1, iter = bayes_iter, 
                  warmup = bayes_warmup, seed = 12345 + j, refresh = 0, silent = 2), silent = TRUE)
                post_samples <- get_valid_group_draws(fit_bayes_interim)
                if (is.null(post_samples)) {
                  simulation_ok <- FALSE
                  break
                }
                interim_PP[k] <- mean(post_samples < 0)
            }
            if (!simulation_ok) {
                return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
            }
            df_final_clean <- df_j[, c("time", "status", "group")]
            fit_bayes_final <- try(update(base_model, newdata = df_final_clean, chains = 4, cores = 1, iter = bayes_iter, 
                warmup = bayes_warmup, seed = 12345 + j, refresh = 0, silent = 2), silent = TRUE)
            post_samples <- get_valid_group_draws(fit_bayes_final)
            if (is.null(post_samples)) {
                return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
            }
            PP_final <- mean(post_samples < 0)
            list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
        }, mc.cores = num_cores_validate)
        saveRDS(stats_list, cache_file)
    }
    invisible(stats_list)
}
generate_calibration_scenario("H0 [survival_colon]","H0")
generate_calibration_scenario("H1_calibration [standard effect]","H1_calibration")
cat("PP_CACHE_GENERATION_COMPLETE: 10000 attempted paths per H0/H1-standard scenario; parameter selection is separate.\n")
