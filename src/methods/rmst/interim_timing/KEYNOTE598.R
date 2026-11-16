if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Independent Monte Carlo validation.
# Method: RMST; case: KEYNOTE-598
# Locked decision parameters: PP_final_threshold = 0.864;
#                             lambda_prime_bayes = 0.210
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_SENS_LOOKS_50_80"
iv_project_id <- "RMST_KEYNOTE598"
iv_method <- "RMST"
iv_case <- "KEYNOTE-598"
iv_candidate_id <- "4076"
iv_locked_pp <- 0.864
iv_locked_lambda <- 0.210
iv_data_seed_base <- 2822001L
iv_model_seed_base <- 2922001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_dgm_file <- file.path(Sys.getenv("V33_CODE_ROOT"), "src", "dgm", "primary", "KEYNOTE598", "data_simulation_final.R")
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

if (!file.exists(iv_dgm_file)) stop("Missing case-specific DGM: ", iv_dgm_file)
source(
  file.path(iv_version_dir, "rmst_iv_a05_a06_helpers.R"),
  local = TRUE
)

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
    "replicate_output=replicate_id,fit_ok,failure_reason,PP_interims,PP_final,early_futility,stop_stage,sample_size_used,go,final_no_go",
    "A05_output=probability,MCSE,Wilson_95pct_Monte_Carlo_interval",
    "A06_output=direct_final_no_go_and_partition_identity_check",
    "failure_handling=reported_both_valid_only_and_all_replicates_failure_as_no_go",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(iv_output_dir, "RUN_METADATA.txt")
)

# 对 KEYNOTE598 H1_observed 生成数据进行自适应设计验证 (RMST 版本)
# Locked candidate 4076: PP = 0.864, lambda = 0.210.

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

rmst_weibull <- function(mu, alpha, tau) {
  scale <- mu / gamma(1 + 1 / alpha)
  x <- (tau / scale)^alpha
  
  res <- mu * pgamma(x, shape = 1 / alpha, rate = 1)
  
  # 检查 RMST 是否在 0 至 tau 范围内。
  stopifnot(all(res >= 0))
  stopifnot(all(res <= tau + 1e-8))
  
  return(res)
}

# 2. 重复次数与结果输出设置
n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
results_file <- if (N.sim_realcase == 10000) {
  "independent_validation_results_10k_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.txt"
} else {
  "independent_validation_results_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.txt"
}
results_file <- file.path(iv_output_dir, basename(results_file))
sink(results_file, split = TRUE)

source(iv_dgm_file, local = TRUE)
trial_settings <- keynote598_trial_settings()
# Pure independent validation: no preliminary single-dataset analysis.
# 4. 决策参数
# Locked before independent validation; no environment override is permitted.
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- trial_settings$N.max
interim_counts <- c(284L, 454L)
interim_fracs <- as.numeric(interim_counts) / as.numeric(N.max)
stopifnot(
  length(interim_counts) == length(interim_fracs),
  all(interim_counts > 0L),
  all(interim_counts < N.max)
)
recruit.int <- trial_settings$recruit.int
tau_fixed <- c(5, 9, 12)
bayes_iter <- 8000L
bayes_warmup <- 4000L
bayes_adapt <- 0.95
bayes_treedepth <- 12L

# 打乱受试者以模拟真实随机入组顺序
# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# 预编译基础贝叶斯模型 (Weibull 家族)
cat("\n正在预编译基础贝叶斯 Weibull 模型 (以极大提速后续分析)...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = weibull(),
  prior   = c(
    prior(normal(0, 1.5), class = "b", coef = "group"),
    prior(normal(log(20), 1.5), class = "Intercept"),
    prior(lognormal(0, 1), class = "shape")
  ),
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
cat("Base model precompilation completed.\n\n")

# 7. 针对 KEYNOTE598 三场景进行自适应设计模拟验证
cat("\n=========================================================================\n")
cat(sprintf("正在对三场景进行 %d 次自适应设计模拟特征验证 (RMST 法)...\n", N.sim_realcase))
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source(iv_dgm_file, local = TRUE)



cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

# N.sim_realcase is parsed at the top of the file

run_validation_scenario <- function(label, scenario) {
  cache_suffix <- if (N.sim_realcase == 10000) "10k_" else ""
  cache_file <- sprintf(
    "validate_stats_%s%s_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.rds",
    cache_suffix, scenario
  )
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
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      # 过滤 PFS <= 0 的点以防止 Weibull 报错
      df_j <- df_j[df_j$progression_t > 0 & df_j$death_t > 0, ]
      if (nrow(df_j) < 10) {
        return(list(
          ok = FALSE,
          failure_reason = "insufficient_generated_rows",
          PP = rep(NA_real_, length(interim_fracs)),
          PP_final = NA_real_
        ))
      }
      
      time_to_first <- pmin(df_j$death_t, df_j$progression_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      original_censor_t <- df_j$censor_t
      
      # 对每个期中分析做自适应截断和模型拟合
      for (k in seq_along(interim_fracs)) {
    n_k <- interim_counts[k]
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0 & time_interim > 0
        
        df_interim <- data.frame(
          time   = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group  = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_bayes_interim <- try(
            update(
              base_model,
              newdata = df_interim,
              chains  = 4,
              cores   = 1,
              iter    = bayes_iter,
              warmup  = bayes_warmup,
              seed    = iv_model_seed_base + j,
              refresh = 0,
              silent  = 2,
              control = list(
                adapt_delta   = bayes_adapt,
                max_treedepth = bayes_treedepth
              )
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_bayes_interim, "try-error")) {
            newdat_sim <- data.frame(group = c(0, 1))
            mu_draws_sim <- try(posterior_linpred(fit_bayes_interim, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
            alpha_draws_sim <- try(posterior_linpred(fit_bayes_interim, dpar = "shape", transform = TRUE), silent = TRUE)
            
            if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
              mu0_draws <- mu_draws_sim[, 1]
              mu1_draws <- mu_draws_sim[, 2]
              alpha_vec <- alpha_draws_sim[, 1]
              
              rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_k)
              rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_k)
              
              interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
            }
          }
        }
      }
      
      # 最终分析 (施加试验启动算起 41 个月的行政删失日历时间截断)
      Time.final <- trial_settings$final_analysis_time
      obs_window_final <- pmax(Time.final - entry_times, 0)
      censor_t_final <- obs_window_final
      
      sub_df <- df_j[, c("time", "status", "group")]
      sub_df$time <- pmin(df_j$time, censor_t_final)
      sub_df$status <- as.integer(df_j$status == 1 & df_j$time <= censor_t_final)
      tau_final <- tau_fixed[length(tau_fixed)]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_bayes_final <- try(
          update(
            base_model,
            newdata = sub_df,
            chains  = 4,
            cores   = 1,
            iter    = bayes_iter,
            warmup  = bayes_warmup,
            seed    = iv_model_seed_base + j,
            refresh = 0,
            silent  = 2,
            control = list(
              adapt_delta   = bayes_adapt,
              max_treedepth = bayes_treedepth
            )
          ),
          silent = TRUE
        )
        if (!inherits(fit_bayes_final, "try-error")) {
          newdat_sim <- data.frame(group = c(0, 1))
          mu_draws_sim <- try(posterior_linpred(fit_bayes_final, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
          alpha_draws_sim <- try(posterior_linpred(fit_bayes_final, dpar = "shape", transform = TRUE), silent = TRUE)
          
          if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
            mu0_draws <- mu_draws_sim[, 1]
            mu1_draws <- mu_draws_sim[, 2]
            alpha_vec <- alpha_draws_sim[, 1]
            
            rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_final)
            rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_final)
            
            PP_final <- mean(rmst1_draws > rmst0_draws)
          }
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else {
          "mcmc_fit_or_posterior_extraction_failed"
        },
        PP = interim_PP,
        PP_final = PP_final
      )
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  c_primes <- PP_final_threshold * (interim_counts / N.max)^lambda_prime_bayes

  iv_a05_a06 <- rmst_iv_write_a05_a06_outputs(
    stats_list = stats_list,
    label = label,
    c_primes = c_primes,
    N.max = N.max,
    interim_fracs = interim_fracs,
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    output_dir = iv_output_dir
  )
  valid_stats <- stats_list[iv_a05_a06$valid_indices]
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的RMST拟合结果。", label))
  }
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(interim_counts[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf(
    "  - 有效模拟次数             : %d/%d\n",
    iv_a05_a06$n_valid,
    iv_a05_a06$n_total
  ))
  cat(sprintf(
    "  - 拟合失败次数             : %d\n",
    iv_a05_a06$n_failed
  ))
  cat("  - A05输出                  : 概率、MCSE、Wilson 95% Monte Carlo区间\n")
  cat("  - A06输出                  : final no-go直接计数及PET+go+final no-go恒等核对\n")
  cat(sprintf("  - 预期早停率 (PET)     : %.2f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.2f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [KEYNOTE598]", "H0")
run_validation_scenario("H1_calibration [standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [pembrolizumab + ipilimumab observed effect]", "H1_observed")

sink()
cat(sprintf(
  "\nInterim-timing sensitivity analysis completed. Results: %s\n",
  results_file
))
