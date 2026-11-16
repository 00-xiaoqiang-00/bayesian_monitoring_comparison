if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Independent Monte Carlo validation.
# Method: Cox; case: IMpower131
# Decision parameters are hard-coded from locked candidate 9610.
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# KM-target rho sensitivity validation; no recalibration or decision-parameter selection.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_SENSITIVITY_RHO_KM_TARGET"
iv_project_id <- "cox_IMpower"
iv_method <- "Cox"
iv_case <- "IMpower131"
iv_candidate_id <- "9610"
iv_locked_pp <- 0.859
iv_locked_lambda <- 0.250
if (!is.finite(iv_locked_pp) || iv_locked_pp <= 0 || iv_locked_pp >= 1) {
  stop("Set LOCKED_PP_FINAL_THRESHOLD to the parameter locked in the supplied calibration settings.")
}
if (!is.finite(iv_locked_lambda) || iv_locked_lambda < 0) {
  stop("Set LOCKED_LAMBDA_PRIME to the parameter locked in the supplied calibration settings.")
}
iv_baseline_data_seed_base <- 1811001L
iv_baseline_model_seed_base <- 1911001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

iv_baseline_rho <- 0.5
iv_allowed_rho <- c(low = 0.3, high = 0.7)
iv_rho_level <- Sys.getenv("SENSITIVITY_RHO_LEVEL", unset = "")
iv_rho_raw <- Sys.getenv("SENSITIVITY_RHO", unset = "")
if (!nzchar(iv_rho_level) || !iv_rho_level %in% names(iv_allowed_rho)) {
  stop("SENSITIVITY_RHO_LEVEL is required and must be exactly 'low' or 'high'.")
}
if (!nzchar(iv_rho_raw)) {
  stop("SENSITIVITY_RHO is required.")
}
iv_sensitivity_rho <- suppressWarnings(as.numeric(iv_rho_raw))
if (length(iv_sensitivity_rho) != 1L || !is.finite(iv_sensitivity_rho) ||
    abs(iv_sensitivity_rho) >= 1) {
  stop("SENSITIVITY_RHO must be one finite numeric value with abs(rho) < 1.")
}
iv_expected_rho <- unname(iv_allowed_rho[[iv_rho_level]])
if (!isTRUE(all.equal(iv_sensitivity_rho, iv_expected_rho, tolerance = 1e-12))) {
  stop(sprintf(
    "SENSITIVITY_RHO mismatch: level '%s' requires rho=%s for %s.",
    iv_rho_level, format(iv_expected_rho, digits = 15), iv_case
  ))
}
iv_rho_safe_tag <- gsub(".", "p", format(
  iv_sensitivity_rho, scientific = FALSE, trim = TRUE, digits = 15
), fixed = TRUE)
iv_rho_safe_tag <- gsub("-", "m", iv_rho_safe_tag, fixed = TRUE)
iv_fixed_interim_schedule <- "272,408,544 (40%,60%,80%)"

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

iv_seed_offset <- 0L
# Seed bases are set by the KM-target helper after project data-generation code is sourced.

iv_output_root <- Sys.getenv("RHO_SENSITIVITY_OUTPUT_ROOT", unset = Sys.getenv("V33_OUTPUT_DIR"))
iv_output_dir <- file.path(
  iv_output_root, "outputs", iv_method, iv_project_id, paste0("rho_", iv_rho_safe_tag), iv_run_tag
)
iv_cache_dir <- file.path(iv_output_dir, "cache")
dir.create(iv_cache_dir, recursive = TRUE, showWarnings = FALSE)

iv_write_completion_marker <- function() {
  if (sink.number(type = "output") != 0L) {
    stop("Refusing to write RUN_COMPLETE.txt while an output sink is active.")
  }
  marker_file <- file.path(iv_output_dir, "RUN_COMPLETE.txt")
  temporary_file <- tempfile(
    pattern = ".RUN_COMPLETE_", tmpdir = iv_output_dir, fileext = ".tmp"
  )
  on.exit({
    if (file.exists(temporary_file)) unlink(temporary_file)
  }, add = TRUE)
  writeLines(
    c(
      paste0("version=", iv_version),
      paste0("project=", iv_project_id),
      paste0("rho_level=", iv_rho_level),
      paste0("sensitivity_rho=", format(iv_sensitivity_rho, digits = 15)),
      paste0("run_tag=", iv_run_tag),
      paste0("completed_at=", format(
        Sys.time(), tz = "Asia/Shanghai", usetz = TRUE
      ))
    ),
    temporary_file
  )
  if (!isTRUE(file.rename(temporary_file, marker_file))) {
    stop("Failed to atomically rename completion marker: ", marker_file)
  }
  invisible(marker_file)
}

iv_seed_from_label <- function(label) {
  label_codes <- utf8ToInt(enc2utf8(as.character(label)))
  label_offset <- if (length(label_codes)) {
    sum(seq_along(label_codes) * label_codes)
  } else {
    0
  }
  as.integer((iv_data_seed_base + label_offset) %% .Machine$integer.max)
}

# RUN_METADATA.txt is written after KM-target parameters and seeds are locked.

# 对 IMpower131 真实重建数据进行自适应设计验证
# 锁定参数: PP_final_threshold = 0.859, lambda_prime_bayes = 0.250

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

# Simulated data are generated by the case DGM.
# 4. 决策参数
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- 680
interim_fracs <- c(0.40, 0.60, 0.80)
recruit.int <- 22 / 680 # 入组速度 (保持与之前校准一致)
bayes_iter <- 2000
bayes_warmup <- 1000 # 每链保留 1000；4 链共 4000 个 post-warmup draws
cache_tag <- "km_target_rho_sensitivity_commonCutoff41_iter2000_warmup1000"

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

# 按照随机化临床试验的实际情况，打乱受试者入组顺序
# 随机打乱分组后用于中期分析。
# Planned entry times used by all newly generated validation replicates.
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

# 7. 针对两个敏感性场景（H0 与 H1 标准场景）进行大规模自适应设计验证
cat("\n=========================================================================\n")
cat("正在对两个敏感性场景进行自适应设计模拟验证 (Cox 法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
iv_data_generator_file <- Sys.getenv("RHO_SENSITIVITY_DATA_GENERATOR", unset = "")
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
# 删失感知参数。第一维必须是复发时间，不是 PFS。

# Rho sensitivity changes only the latent log-time covariance off-diagonal.
iv_rebuild_sigma_for_rho <- function(Sigma_baseline, rho, label) {
  if (!is.matrix(Sigma_baseline) || !identical(dim(Sigma_baseline), c(2L, 2L)) ||
      any(!is.finite(Sigma_baseline)) || any(diag(Sigma_baseline) <= 0)) {
    stop(label, " baseline covariance must be a finite 2x2 matrix with positive diagonal.")
  }
  Sigma <- Sigma_baseline
  covariance <- rho * sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  Sigma[1L, 2L] <- covariance
  Sigma[2L, 1L] <- covariance
  if (!identical(diag(Sigma), diag(Sigma_baseline)) ||
      min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " rho covariance reconstruction failed.")
  }
  Sigma
}

iv_assert_baseline_sigma <- function(Sigma_baseline, label) {
  expected <- iv_baseline_rho * sqrt(prod(diag(Sigma_baseline)))
  if (abs(Sigma_baseline[1L, 2L] - expected) > 5e-5 ||
      abs(Sigma_baseline[2L, 1L] - expected) > 5e-5) {
    stop(label, " does not reconstruct the locked baseline rho.")
  }
  invisible(TRUE)
}

mu.ctrl.baseline <- c(log(6.1202), log(13.4685))   # 对照组：复发时间 / 死亡时间
mu.ctrl <- mu.ctrl.baseline
stopifnot(identical(mu.ctrl, mu.ctrl.baseline))
Sigma.ctrl.baseline <- matrix(c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2)
iv_assert_baseline_sigma(Sigma.ctrl.baseline, "IMpower control")
Sigma.ctrl <- iv_rebuild_sigma_for_rho(
  Sigma.ctrl.baseline, iv_sensitivity_rho, "IMpower control"
)

mu.trt.real <- c(log(7.5610), log(14.3431))   # 试验组：复发时间 / 死亡时间
Sigma.trt.real <- matrix(c(0.8207, 0.5298, 0.5298, 1.3679), nrow = 2)

# 标准疗效 H1 场景（复发时间=8.00月，死亡时间=18.00月）
mu.trt.H1_standard.baseline <- c(log(8.00), log(18.00))
mu.trt.H1_standard <- mu.trt.H1_standard.baseline
stopifnot(identical(mu.trt.H1_standard, mu.trt.H1_standard.baseline))
Sigma.trt.H1_standard.baseline <- matrix(
  c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2
)
iv_assert_baseline_sigma(Sigma.trt.H1_standard.baseline, "IMpower H1-standard")
Sigma.trt.H1_standard <- iv_rebuild_sigma_for_rho(
  Sigma.trt.H1_standard.baseline, iv_sensitivity_rho, "IMpower H1-standard"
)


# Load case- and scenario-specific KM-target parameters and seed rules.
iv_rho_helper <- Sys.getenv("RHO_SENSITIVITY_HELPER", unset = "")
if (!nzchar(iv_rho_helper)) stop("RHO_SENSITIVITY_HELPER is required.")
source(iv_rho_helper, local = TRUE)
iv_master_seed <- as.integer(Sys.getenv("RHO_SENSITIVITY_MASTER_SEED", unset = "8252026"))
if (!is.finite(iv_master_seed)) stop("Invalid RHO_SENSITIVITY_MASTER_SEED.")
iv_data_seed_base <- rho_data_seed_base(iv_master_seed, iv_case, iv_rho_level)
iv_model_seed_base <- rho_model_seed_base(iv_master_seed, iv_case, iv_method, iv_rho_level)
if (iv_smoke_test) {
  iv_data_seed_base <- as.integer((iv_data_seed_base + 500000000L) %% (.Machine$integer.max - 1))
  iv_model_seed_base <- as.integer((iv_model_seed_base + 700000000L) %% (.Machine$integer.max - 1))
}
iv_formal_data_seed_base <- rho_data_seed_base(iv_master_seed, iv_case, iv_rho_level)
iv_formal_model_seed_base <- rho_model_seed_base(iv_master_seed, iv_case, iv_method, iv_rho_level)
iv_sensitivity_formal_data_seed_base <- iv_formal_data_seed_base
iv_sensitivity_formal_model_seed_base <- iv_formal_model_seed_base
iv_seed_from_label <- function(label) {
  scenario_tag <- if (grepl("H0", label, fixed = TRUE)) "H0" else if (
    grepl("H1_observed|Realcase", label)
  ) "H1_observed" else "H1_standard"
  rho_data_seed(iv_master_seed, iv_case, scenario_tag, iv_rho_level)
}
iv_parameter_table_sha256 <- Sys.getenv("RHO_SENSITIVITY_PARAMETER_SHA256", unset = "not_supplied")
iv_parameter_table_path <- Sys.getenv("RHO_PARAMETER_FILE")
iv_write_km_metadata <- function() {
  writeLines(c(
    paste0("version=", iv_version), paste0("project=", iv_project_id),
    paste0("method=", iv_method), paste0("case=", iv_case),
    paste0("candidate_id=", iv_candidate_id), paste0("run_tag=", iv_run_tag),
    "sensitivity_only=true", "recalibration=false",
    "sensitivity_parameter=KM_target_refitted_latent_DGM_rho",
    "parameter_source=locked_KM_target_parameter_table",
    paste0("parameter_table=", iv_parameter_table_path),
    paste0("parameter_table_sha256=", iv_parameter_table_sha256),
    paste0("sensitivity_rho=", format(iv_sensitivity_rho, digits = 15)),
    paste0("baseline_rho=", format(iv_baseline_rho, digits = 15)),
    paste0("rho_level=", iv_rho_level),
    paste0("fixed_interim_schedule=", iv_fixed_interim_schedule),
    "scenarios=H0,H1_standard,H1_observed",
    "baseline_decision_parameters_locked=true",
    "realcase_or_H1_observed_not_used_for_parameter_selection=true",
    "cache_isolation=KM-target rho-specific output directory",
    paste0("locked_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("locked_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    "locked_RMUT_weight=0.5",
    paste0("master_seed=", iv_master_seed),
    paste0("formal_data_seed_base=", iv_formal_data_seed_base),
    paste0("formal_model_seed_base=", iv_formal_model_seed_base),
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    "seed_policy=data stream is shared across methods for case/scenario/rho; model stream is method-specific; smoke uses explicit offsets",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ), file.path(iv_output_dir, "RUN_METADATA.txt"))
}

iv_ctrl_km <- rho_make_impower_parameters("H0", "control", iv_rho_level)
iv_ctrl_h1_standard_km <- rho_make_impower_parameters("H1-standard", "control", iv_rho_level)
iv_ctrl_h1_observed_km <- rho_make_impower_parameters("H1-observed", "control", iv_rho_level)
iv_h1_standard_km <- rho_make_impower_parameters("H1-standard", "treatment", iv_rho_level)
iv_h1_observed_km <- rho_make_impower_parameters("H1-observed", "treatment", iv_rho_level)
mu.ctrl <- iv_ctrl_km$mu; Sigma.ctrl <- iv_ctrl_km$Sigma
mu.ctrl.H1_standard <- iv_ctrl_h1_standard_km$mu; Sigma.ctrl.H1_standard <- iv_ctrl_h1_standard_km$Sigma
mu.ctrl.H1_observed <- iv_ctrl_h1_observed_km$mu; Sigma.ctrl.H1_observed <- iv_ctrl_h1_observed_km$Sigma
mu.trt.H1_standard <- iv_h1_standard_km$mu; Sigma.trt.H1_standard <- iv_h1_standard_km$Sigma
mu.trt.H1 <- iv_h1_standard_km$mu; Sigma.trt.H1 <- iv_h1_standard_km$Sigma
mu.trt.real <- iv_h1_observed_km$mu; Sigma.trt.real <- iv_h1_observed_km$Sigma
mu.trt.H1_sens2 <- mu.trt.H1_standard; Sigma.trt.H1_sens2 <- Sigma.trt.H1_standard
iv_write_km_metadata()

run_validation_scenario <- function(label, mu_trt, Sigma_trt, mu_ctrl, Sigma_ctrl) {
  scenario_tag <- ifelse(grepl("H0", label), "H0", ifelse(grepl("H1_observed", label), "H1_observed", "H1_standard"))
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
      mu.trt = mu_trt, Sigma.trt = Sigma_trt, 
      mu.ctrl = mu_ctrl, Sigma.ctrl = Sigma_ctrl, cens_upper = 1e6
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
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
      
      final_admin_followup <- pmax(41 - entry_times, 0)
      final_censor_t <- pmin(original_censor_t, final_admin_followup)
      final_time <- pmin(time_to_first, final_censor_t)
      final_status <- as.integer(time_to_first <= final_censor_t)
      final_valid_rows <- final_censor_t > 0
      df_final_clean <- data.frame(
        time = final_time[final_valid_rows],
        status = final_status[final_valid_rows],
        group = df_j$group[final_valid_rows]
      )
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
  cat(sprintf("  - 预期早停率 (PET)     : %.2f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.2f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 (原假设/药物无效) [校准用场景]", mu.ctrl, Sigma.ctrl, mu.ctrl, Sigma.ctrl)
run_validation_scenario("H1_standard [KM-target standard effect]", mu.trt.H1_standard, Sigma.trt.H1_standard, mu.ctrl.H1_standard, Sigma.ctrl.H1_standard)
run_validation_scenario("H1_observed [KM-target observed effect]", mu.trt.real, Sigma.trt.real, mu.ctrl.H1_observed, Sigma.ctrl.H1_observed)

sink()

iv_write_completion_marker()
