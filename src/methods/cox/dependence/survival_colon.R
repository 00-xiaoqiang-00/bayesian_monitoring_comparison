if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Independent Monte Carlo validation.
# Method: Cox; case: survival_colon
# Locked decision parameters: PP_final_threshold = 0.851;
#                             lambda_prime_bayes = 0.285
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# KM-target rho sensitivity validation; no recalibration or decision-parameter selection.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_SENSITIVITY_RHO_KM_TARGET"
iv_project_id <- "cox_survival_colon"
iv_method <- "Cox"
iv_case <- "survival_colon"
iv_candidate_id <- "10939"
iv_locked_pp <- 0.851
iv_locked_lambda <- 0.285
iv_baseline_data_seed_base <- 1814001L
iv_baseline_model_seed_base <- 1914001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

iv_baseline_rho <- 0.9082540001403667
iv_allowed_rho <- c(low = 0.7, high = 0.98)
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
iv_fixed_interim_schedule <- "248,372,496 (40%,60%,80%)"

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

# 对 survival_colon H1_observed 生成数据进行自适应设计验证
# Locked candidate 10939: PP = 0.851, lambda = 0.285.

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

iv_data_generator_file <- Sys.getenv("RHO_SENSITIVITY_DATA_GENERATOR", unset = "")
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)
# Generate and analyse the scenario replicates.
trial_settings <- survival_colon_trial_settings()
# 4. 决策参数
# Locked before independent validation; no environment override is permitted.
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- trial_settings$N.max
interim_fracs <- trial_settings$interim_fracs
recruit.int <- trial_settings$recruit.int
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

# Run independent validation with the case data generator.
# 7. 针对 survival_colon 两个敏感性场景进行大规模自适应设计验证
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

# Rebuild only control/H1-calibration covariance off-diagonals at candidate rho.
iv_rebuild_colon_sigma <- function(Sigma_baseline, rho, label) {
  if (!is.matrix(Sigma_baseline) || !identical(dim(Sigma_baseline), c(2L, 2L)) ||
      any(!is.finite(Sigma_baseline)) || any(diag(Sigma_baseline) <= 0)) {
    stop(label, " baseline covariance is invalid.")
  }
  baseline_reconstructed_rho <- Sigma_baseline[1L, 2L] /
    sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  if (!isTRUE(all.equal(baseline_reconstructed_rho, iv_baseline_rho,
                        tolerance = 1e-12)) ||
      !isTRUE(all.equal(Sigma_baseline[2L, 1L], Sigma_baseline[1L, 2L],
                        tolerance = 1e-12))) {
    stop(label, " does not reconstruct the declared locked baseline rho.")
  }
  Sigma <- Sigma_baseline
  covariance <- rho * sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  Sigma[1L, 2L] <- covariance
  Sigma[2L, 1L] <- covariance
  if (!identical(diag(Sigma), diag(Sigma_baseline)) ||
      min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " candidate covariance is not positive-definite.")
  }
  Sigma
}

iv_baseline_parameters <- survival_colon_generator_parameters()
iv_sensitivity_parameters <- iv_baseline_parameters
for (iv_arm_name in c("control", "H1_calibration")) {
  iv_baseline_arm <- iv_baseline_parameters[[iv_arm_name]]
  iv_sensitivity_parameters[[iv_arm_name]]$Sigma <- iv_rebuild_colon_sigma(
    iv_baseline_arm$Sigma, iv_sensitivity_rho, iv_arm_name
  )
  if (!identical(iv_sensitivity_parameters[[iv_arm_name]]$mu,
                 iv_baseline_arm$mu) ||
      !identical(diag(iv_sensitivity_parameters[[iv_arm_name]]$Sigma),
                 diag(iv_baseline_arm$Sigma))) {
    stop(iv_arm_name, " survival_colon marginal parameter changed.")
  }
  survival_colon_validate_arm(iv_sensitivity_parameters[[iv_arm_name]], iv_arm_name)
}
if (!identical(iv_sensitivity_parameters$H1_observed,
               iv_baseline_parameters$H1_observed) ||
    !identical(iv_sensitivity_parameters$targets,
               iv_baseline_parameters$targets)) {
  stop("Non-target survival_colon parameters changed.")
}


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

iv_sensitivity_parameters_by_scenario <- rho_make_generator_parameters_by_scenario(iv_case, iv_rho_level)
iv_sensitivity_parameters <- iv_sensitivity_parameters_by_scenario$H0
if (iv_method == "RMUT" && nzchar(Sys.getenv("RMUT_WEIGHT", unset = ""))) {
  stop("RMUT_WEIGHT override is forbidden in KM-target rho sensitivity validation.")
}
iv_write_km_metadata()

run_validation_scenario <- function(label, scenario) {
  scenario_tag <- if (identical(scenario, "H1_calibration")) "H1_standard" else scenario
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
      cens_upper = trial_settings$cens_upper,
      scenario = scenario,
      parameters = iv_sensitivity_parameters_by_scenario[[scenario]]
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      final_censor_t <- pmin(
        df_j$censor_t,
        survival_colon_admin_followup(
          entry_times, trial_settings$final_analysis_time
        )
      )
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
  cat(sprintf("  - 预期早停率 (PET)     : %.2f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.2f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [survival_colon]", "H0")
run_validation_scenario("H1_standard [KM-target standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [KM-target observed effect]", "H1_observed")

sink()

iv_write_completion_marker()
