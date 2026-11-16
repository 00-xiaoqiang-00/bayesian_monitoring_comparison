if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Independent Monte Carlo validation.
# Method: Cox; case: OAK
# Locked decision parameters: PP_final_threshold = 0.850;
#                             lambda_prime_bayes = 0.340
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_STANDARD"
iv_project_id <- "cox_OAK"
iv_method <- "Cox"
iv_case <- "OAK"
iv_candidate_id <- "13039"
iv_locked_pp <- 0.850
iv_locked_lambda <- 0.340
iv_data_seed_base <- 1813001L
iv_model_seed_base <- 1913001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

iv_smoke_test <- identical(Sys.getenv("IV_SMOKE_TEST", unset = "0"), "1")
if (iv_smoke_test) {
  iv_smoke_n <- suppressWarnings(as.integer(Sys.getenv("N_SIM_REALCASE", unset = "")))
  if (!is.finite(iv_smoke_n) || iv_smoke_n < 1L || iv_smoke_n > 20L) {
    stop("Smoke mode requires N_SIM_REALCASE between 1 and 20.")
  }
  iv_run_tag <- sprintf("smoke_n%d", iv_smoke_n)
} else {
  Sys.setenv(N_SIM_REALCASE = "10000")
  iv_run_tag <- "formal_10000"
}

iv_formal_data_seed_base <- iv_data_seed_base
iv_formal_model_seed_base <- iv_model_seed_base
if (iv_smoke_test) {
  iv_data_seed_base <- iv_data_seed_base + 500000000L
  iv_model_seed_base <- iv_model_seed_base + 500000000L
}

iv_output_dir <- file.path(Sys.getenv("V33_OUTPUT_DIR"), iv_run_tag)
iv_cache_dir <- file.path(iv_output_dir, "cache")
dir.create(iv_cache_dir, recursive = TRUE, showWarnings = FALSE)

iv_seed_from_label <- function(label) {
  label_codes <- utf8ToInt(enc2utf8(as.character(label)))
  label_offset <- if (length(label_codes)) {
    sum(seq_along(label_codes) * label_codes)
  } else {
    0
  }
  as.integer((iv_data_seed_base + label_offset) %% .Machine$integer.max)
}

writeLines(
  c(
    paste0("version=", iv_version),
    paste0("project=", iv_project_id),
    paste0("method=", iv_method),
    paste0("case=", iv_case),
    paste0("candidate_id=", iv_candidate_id),
    paste0("run_tag=", iv_run_tag),
    paste0("locked_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("locked_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    paste0("formal_data_seed_base=", iv_formal_data_seed_base),
    paste0("formal_model_seed_base=", iv_formal_model_seed_base),
    "seed_policy=smoke test uses +500000000 offset; formal streams remain untouched",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(iv_output_dir, "RUN_METADATA.txt")
)

# independent_validation 案例的独立模拟验证。
# Decision thresholds use this design's configured parameters.

# 1. 加载依赖并配置编译环境
# Optional compiler settings may be supplied through R_MAKEVARS_USER.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- file.path(dirname(dirname(R.home())), "bin")
Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(survival)
library(brms)
library(posterior)

# 2. 分析配置
sink(file.path(iv_output_dir, "independent_validation_results.txt"), split = TRUE)

source("data_simulation_final.R", local = TRUE)
# Pure independent validation: no preliminary single-dataset analysis.
# 4. 配置自适应设计参数
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- 850
interim_fracs <- c(0.40, 0.60, 0.80)
recruit.int <- 8.6 / 850
bayes_iter <- 2000
bayes_warmup <- 1000 # 每链保留 1000；4 链共 4000 个 post-warmup draws
cache_tag <- "OAK_mix_v4_H0_H1calibration_H1observed_cox_iter2000_warmup1000"

get_valid_group_draws <- function(fit) {
  if (inherits(fit, "try-error")) return(NULL)

  draws <- try(posterior::as_draws_array(fit, variable = "b_group"), silent = TRUE)
  if (inherits(draws, "try-error")) return(NULL)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NULL)
  }

  posterior::as_draws_matrix(draws)[, "b_group"]
}

# Planned entry times follow the configured recruitment interval.
entry_times <- recruit.int * (0:(N.max - 1))

# 预编译基础贝叶斯模型
cat("\n正在预编译基础贝叶斯模型 (以极大提速后续分析)...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = brmsfamily("cox"),
  prior   = prior(normal(0, 1.5), class = "b", coef = "group"),
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
cat("Base model precompilation completed.\n\n")

# Generate the scenario replicates.
# 7. 针对 OAK 三场景进行大规模自适应设计验证
cat("\n=========================================================================\n")
cat("正在对三场景进行自适应设计模拟验证 (Cox 法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source("data_simulation_final.R", local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
# 删失感知参数。第一维必须是复发时间，不是 PFS。
run_validation_scenario <- function(label, scenario) {
  scenario_tag <- scenario
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = 109.3,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      final_censor_t <- pmin(df_j$censor_t, oak_admin_followup(entry_times, 27.9))
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
        df_interim <- data.frame(
          time = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        fit_bayes_interim <- try(
          update(base_model, newdata = df_interim, chains = 4, cores = 1,
                 iter = bayes_iter, warmup = bayes_warmup,
                 seed = iv_model_seed_base + j, refresh = 0, silent = 2),
          silent = TRUE
        )
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
      fit_bayes_final <- try(
        update(base_model, newdata = df_final_clean, chains = 4, cores = 1,
               iter = bayes_iter, warmup = bayes_warmup,
               seed = iv_model_seed_base + j, refresh = 0, silent = 2),
        silent = TRUE
      )
      post_samples <- get_valid_group_draws(fit_bayes_final)
      if (is.null(post_samples)) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      PP_final <- mean(post_samples < 0)
      
      list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  stats_list <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list)
  if (length(stats_list) == 0) stop("没有通过 MCMC 诊断的有效验证模拟。")
  cat(sprintf("  - 通过 MCMC 诊断的有效模拟: %d/%d\n",
              length(stats_list), N.sim_realcase))

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(stats_list, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(stats_list, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(stats_list, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [OAK calibration]", "H0")
run_validation_scenario("H1_calibration [OAK target effect]", "H1_calibration")
run_validation_scenario("H1_observed [OAK published effect]", "H1_observed")

sink()
