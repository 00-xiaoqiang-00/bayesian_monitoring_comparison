if (!nzchar(Sys.getenv("V33_CODE_ROOT")) || !nzchar(Sys.getenv("V33_OUTPUT_DIR"))) stop("Run the registered workflow with an external --output-dir.")
# Independent Monte Carlo validation.
# Method: RMUT; case: IMpower131
# Current defaults are recorded in config/rmut_parameters.tsv; a verified RMUT_LOCKED_PARAMETER_MANIFEST may override them.
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# KM-target rho sensitivity validation; no recalibration or decision-parameter selection.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_RHO_KM_TARGET_ROWPAIR_20261010"
iv_project_id <- "RMUT_IMpower131"
iv_method <- "RMUT"
iv_case <- "IMpower131"
iv_candidate_id <- "1310"
iv_locked_pp <- .883
iv_locked_lambda <- .015
if (!is.finite(iv_locked_pp) || iv_locked_pp <= 0 || iv_locked_pp >= 1) {
  stop("Set LOCKED_PP_FINAL_THRESHOLD to the parameter locked in the supplied calibration settings.")
}
if (!is.finite(iv_locked_lambda) || iv_locked_lambda < 0) {
  stop("Set LOCKED_LAMBDA_PRIME to the parameter locked in the supplied calibration settings.")
}
iv_baseline_data_seed_base <- 1831001L
iv_baseline_model_seed_base <- 1931001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- iv_version_dir
iv_dgm_file <- file.path(Sys.getenv("V33_CODE_ROOT"), "src", "dgm", "primary", "IMpower", "data_simulation_final.R")
source(file.path(iv_version_dir, "..", "..", "common", "rmut_pairing.R"), local = TRUE)
source(file.path(iv_version_dir, "..", "..", "common", "rmut_sensitivity_runtime.R"), local = TRUE)
iv_case_code <- 1L
iv_parameters <- rmut_apply_parameter_manifest(iv_case, list(candidate_id = iv_candidate_id, c = iv_locked_pp, gamma = iv_locked_lambda))
iv_candidate_id <- iv_parameters$candidate_id
iv_locked_pp <- iv_parameters$c
iv_locked_lambda <- iv_parameters$gamma
iv_parameter_source <- iv_parameters$source
Sys.setenv(RHO_SENSITIVITY_VERSION_DIR = normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE))
rmut_preflight(iv_version_dir, is_rho = TRUE)
if ("--preflight-only" %in% commandArgs(trailingOnly = TRUE)) {
  cat("PREFLIGHT_OK ", iv_case, " ", iv_version, " parameter_source=", iv_parameter_source, "\n", sep = "")
  quit(save = "no", status = 0L)
}
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
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
  iv_output_root, "outputs", iv_version, iv_method, iv_project_id, paste0("rho_", iv_rho_safe_tag), iv_run_tag
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
    paste0("pairing_implementation=", rmut_implementation_version),
    paste0("parameter_source=", iv_parameter_source),
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

# IMpower 案例的独立模拟验证。
# Decision thresholds are supplied by the parameter manifest.

# 1. 加载依赖并配置编译环境
# Use the host compiler configuration.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- file.path(dirname(dirname(R.home())), "bin")
Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(survival)
library(brms)
library(posterior)

# 2. 分析配置
sink(file.path(iv_output_dir, "independent_validation_results.txt"), split = TRUE, append = TRUE)

# Simulated data are generated by the case DGM.
# 4. 配置自适应设计参数；使用本文件锁定值
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- 680
interim_fracs <- c(0.40, 0.60, 0.80)
recruit.int <- 22 / 680
bayes_iter <- 2000
bayes_warmup <- 1000
bayes_chains <- 4L
bayes_adapt <- 0.95
bayes_treedepth <- 12L
bayes_retry_iter <- 4000L
bayes_retry_warmup <- 2000L
bayes_retry_adapt <- 0.99
bayes_retry_treedepth <- 15L
if (nzchar(Sys.getenv("RMUT_WEIGHT", unset = ""))) {
  stop("RMUT_WEIGHT overrides are forbidden: rho sensitivity locks w_R=0.5.")
}
w_R          <- 0.5
tau_fixed    <- c(6, 9, 12, 18)
stopifnot(length(tau_fixed) == length(interim_fracs) + 1L)
stopifnot(length(w_R) == 1L, is.finite(w_R), w_R >= 0, w_R <= 1, all(is.finite(tau_fixed)), all(tau_fixed > 0), all(diff(tau_fixed) > 0))

prior_death <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(log(12), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)
prior_composite <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(log(6), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)

# 按照随机化临床试验的实际情况，打乱受试者入组顺序
# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# RMUT 后验计算函数
source(file.path(iv_version_dir, "..", "..", "common", "rmut_pairing.R"), local = TRUE)
# 仅使用通过收敛和有效样本量筛查的后验抽样。
extract_valid_aft_draws <- function(fit) {
  if (inherits(fit, "try-error")) return(NULL)

  variables <- c("b_Intercept", "b_group", "sigma")
  draws <- try(posterior::as_draws_array(fit, variable = variables), silent = TRUE)
  if (inherits(draws, "try-error")) return(NULL)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NULL)
  }

  out <- as.data.frame(fit)[, variables, drop = FALSE]
  attr(out, "rmut_posterior_diagnostics") <- as.data.frame(diagnostics)
  out
}

source(
  file.path(iv_version_dir, "..", "..", "common", "rmut_stable_fit_validation.R"),
  local = TRUE
)

# 分别预编译两个组成模型，以使用终点特异的截距先验。
cat("\n正在预编译死亡和首次复发或死亡 AFT 模型...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model_death <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = lognormal(),
  prior   = prior_death,
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
base_model_composite <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = lognormal(),
  prior   = prior_composite,
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 22345,
  refresh = 0,
  silent  = 2
)
cat("Component-model precompilation completed.\n\n")

# 7. 针对两个敏感性场景（H0 与 H1 标准场景）进行独立自适应设计验证
cat("\n=========================================================================\n")
cat("准备对两个敏感性场景进行自适应设计模拟特征验证 (三状态RMUT法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
iv_data_generator_file <- iv_dgm_file
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
analysis_version <- "IV10000_RHO_KM_TARGET_ROWPAIR_20261010"
cache_signature <- sprintf(
  "%s_w%.2f_tau%s",
  analysis_version, w_R, paste(tau_fixed, collapse = "-")
)

# 删失感知参数（通过 survreg 对数正态 AFT 模型估计，正确处理右删失）

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

mu.ctrl.baseline <- c(log(6.1202), log(13.4685))   # 对照组复发/死亡对数位置参数
mu.ctrl <- mu.ctrl.baseline
stopifnot(identical(mu.ctrl, mu.ctrl.baseline))
Sigma.ctrl.baseline <- matrix(c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2)
iv_assert_baseline_sigma(Sigma.ctrl.baseline, "IMpower control")
Sigma.ctrl <- iv_rebuild_sigma_for_rho(
  Sigma.ctrl.baseline, iv_sensitivity_rho, "IMpower control"
)

mu.trt.real <- c(log(7.5610), log(14.3431))   # 试验组复发/死亡对数位置参数
Sigma.trt.real <- matrix(c(0.8207, 0.5298, 0.5298, 1.3679), nrow = 2)

# 标准H1场景（中位复发时间8.00月，死亡时间18.00月）
mu.trt.H1.baseline <- c(log(8.00), log(18.00))
mu.trt.H1 <- mu.trt.H1.baseline
stopifnot(identical(mu.trt.H1, mu.trt.H1.baseline))

Sigma.trt.H1.baseline <- Sigma.ctrl.baseline
iv_assert_baseline_sigma(Sigma.trt.H1.baseline, "IMpower H1-standard")
Sigma.trt.H1 <- iv_rebuild_sigma_for_rho(
  Sigma.trt.H1.baseline, iv_sensitivity_rho, "IMpower H1-standard"
)


# Load case- and scenario-specific KM-target parameters and seed rules.
iv_rho_helper <- file.path(Sys.getenv("V33_CODE_ROOT"), "src", "common", "rho_parameters.R")
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
    paste0("version=", iv_version),
    paste0("pairing_implementation=", rmut_implementation_version),
    paste0("parameter_source=", iv_parameter_source), paste0("project=", iv_project_id),
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
  scenario_id <- if (grepl("H0", label)) "H0" else if (grepl("H1_observed", label)) "H1_observed" else "H1_standard"
  expected_cache_context <- rmut_pp_cache_context(label)
  iv_scenario_code <- rmut_scenario_code(scenario_id)
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_id, cache_signature)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (rmut_cached_stats_valid(stats_list_temp, expected_cache_context, N.sim_realcase, length(interim_fracs))) {
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
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      fit_pair_audit <- vector("list", length(tau_fixed))
      original_censor_t <- df_j$censor_t
      
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(N.max * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        death_event_t <- df_j$death_t[1:n_k]
        death_time    <- pmin(death_event_t, censor_t_interim)
        death_status  <- as.integer(death_event_t <= censor_t_interim)
        
        composite_event_t <- pmin(
          df_j$recurrence_t[1:n_k],
          df_j$death_t[1:n_k]
        )
        composite_time   <- pmin(composite_event_t, censor_t_interim)
        composite_status <- as.integer(composite_event_t <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0
        
        df_death_interim <- data.frame(
          time = death_time[valid_rows],
          status = death_status[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        
        df_composite_interim <- data.frame(
          time = composite_time[valid_rows],
          status = composite_status[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_death_interim <- rmut_fit_with_stable_retry(
            base_model_death,
            df_death_interim,
            iv_model_seed_base + j
          )
          fit_composite_interim <- rmut_fit_with_stable_retry(
            base_model_composite,
            df_composite_interim,
            iv_model_seed_base + j + 100000L
          )
          
          fit_pair_audit[[k]] <- rmut_stage_audit(fit_death_interim, fit_composite_interim,
            rmut_pair_seed(j, k, iv_scenario_code, iv_case_code), tau_k)
          if (!inherits(fit_death_interim, "try-error") && !inherits(fit_composite_interim, "try-error")) {
            draws_d <- extract_valid_aft_draws(fit_death_interim)
            draws_c <- extract_valid_aft_draws(fit_composite_interim)
            
            if (!is.null(draws_d) && !is.null(draws_c) &&
                all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_d)) &&
                all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_c))) {
              
              rmut_ctrl <- compute_rmut_samples(
                mu_D    = draws_d$b_Intercept,
                sigma_D = draws_d$sigma,
                mu_C    = draws_c$b_Intercept,
                sigma_C = draws_c$sigma,
                tau     = tau_k,
                w_R     = w_R,
    pair_seed = rmut_pair_seed(j, k, iv_scenario_code, iv_case_code)
  )
              
              rmut_trt <- compute_rmut_samples(
                mu_D    = draws_d$b_Intercept + draws_d$b_group,
                sigma_D = draws_d$sigma,
                mu_C    = draws_c$b_Intercept + draws_c$b_group,
                sigma_C = draws_c$sigma,
                tau     = tau_k,
                w_R     = w_R,
    pair_seed = rmut_pair_seed(j, k, iv_scenario_code, iv_case_code)
  )
              
              interim_PP[k] <- mean((rmut_trt - rmut_ctrl) > 0)
            }
          }
        }
      }
      
      final_admin_followup <- pmax(41 - entry_times, 0)
      final_censor_t <- pmin(original_censor_t, final_admin_followup)

      death_event_t <- df_j$death_t
      death_time    <- pmin(death_event_t, final_censor_t)
      death_status  <- as.integer(death_event_t <= final_censor_t)
      
      composite_event_t <- pmin(df_j$recurrence_t, df_j$death_t)
      composite_time   <- pmin(composite_event_t, final_censor_t)
      composite_status <- as.integer(composite_event_t <= final_censor_t)
      
      df_death_final <- data.frame(
        time = death_time,
        status = death_status,
        group = df_j$group
      )
      
      df_composite_final <- data.frame(
        time = composite_time,
        status = composite_status,
        group = df_j$group
      )
      
      tau_final <- tau_fixed[4]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_death_final <- rmut_fit_with_stable_retry(
          base_model_death,
          df_death_final,
          iv_model_seed_base + j
        )
        fit_composite_final <- rmut_fit_with_stable_retry(
          base_model_composite,
          df_composite_final,
          iv_model_seed_base + j + 100000L
        )
        
        fit_pair_audit[[length(tau_fixed)]] <- rmut_stage_audit(fit_death_final, fit_composite_final,
          rmut_pair_seed(j, length(tau_fixed), iv_scenario_code, iv_case_code), tau_final)
        if (!inherits(fit_death_final, "try-error") && !inherits(fit_composite_final, "try-error")) {
          draws_d <- extract_valid_aft_draws(fit_death_final)
          draws_c <- extract_valid_aft_draws(fit_composite_final)
          
          if (!is.null(draws_d) && !is.null(draws_c) &&
              all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_d)) &&
              all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_c))) {
            
            rmut_ctrl <- compute_rmut_samples(
              mu_D    = draws_d$b_Intercept,
              sigma_D = draws_d$sigma,
              mu_C    = draws_c$b_Intercept,
              sigma_C = draws_c$sigma,
              tau     = tau_final,
              w_R     = w_R,
    pair_seed = rmut_pair_seed(j, length(tau_fixed), iv_scenario_code, iv_case_code)
  )
            
            rmut_trt <- compute_rmut_samples(
              mu_D    = draws_d$b_Intercept + draws_d$b_group,
              sigma_D = draws_d$sigma,
              mu_C    = draws_c$b_Intercept + draws_c$b_group,
              sigma_C = draws_c$sigma,
              tau     = tau_final,
              w_R     = w_R,
    pair_seed = rmut_pair_seed(j, length(tau_fixed), iv_scenario_code, iv_case_code)
  )
            
            PP_final <- mean((rmut_trt - rmut_ctrl) > 0)
          }
        }
      }
      
      list(replicate_id = as.integer(j), PP = interim_PP, PP_final = PP_final, fit_pair_audit = fit_pair_audit)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    attr(stats_list, "rmut_pp_cache_context") <- expected_cache_context
    if (!rmut_cached_stats_valid(stats_list, expected_cache_context, N.sim_realcase, length(interim_fracs))) {
      stop("Refusing to save incomplete or malformed corrected PP cache.")
    }
    saveRDS(stats_list, cache_file)
  }
  
  valid_stats <- Filter(function(stats) {
    is.list(stats) && length(stats$PP) == length(interim_fracs) &&
      all(is.finite(stats$PP)) && is.finite(stats$PP_final)
  }, stats_list)
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的MCMC模拟结果。", label))
  }
  if (length(valid_stats) < length(stats_list)) {
    cat(sprintf("警告：%s 排除了 %d 次MCMC失败的模拟。\n",
                label, length(stats_list) - length(valid_stats)))
  }

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes)
    !stop_val && stats$PP_final > PP_final_threshold
  }))
  
  formal_summary <- rmut_summarize_stats(stats_list, label, scenario_id, c_primes,
    floor(N.max * interim_fracs))
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 有效模拟次数             : %d/%d\n", length(valid_stats), length(stats_list)))
  cat(sprintf("  - 预期早停率 (PET)     : %.2f%%\n", formal_summary$PET * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", formal_summary$ESS, (formal_summary$ESS / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.2f%%\n", formal_summary$go_rate * 100))
  cat("=========================================================================\n")
}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 (原假设/药物无效) [校准用场景]", mu.ctrl, Sigma.ctrl, mu.ctrl, Sigma.ctrl)
run_validation_scenario("H1_standard [KM-target standard effect]", mu.trt.H1_standard, Sigma.trt.H1_standard, mu.ctrl.H1_standard, Sigma.ctrl.H1_standard)
run_validation_scenario("H1_observed [KM-target observed effect]", mu.trt.real, Sigma.trt.real, mu.ctrl.H1_observed, Sigma.ctrl.H1_observed)

sink()

rmut_write_completion()

