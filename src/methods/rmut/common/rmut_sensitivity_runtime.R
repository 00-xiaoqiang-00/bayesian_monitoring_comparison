rmut_shared_paths <- function(directory) {
  directory <- normalizePath(directory, mustWork = TRUE)
  common <- normalizePath(file.path(directory, "..", "..", "common"), mustWork = TRUE)
  case <- basename(directory)
  scope <- basename(dirname(directory))
  stable <- if (scope == "calibration") "rmut_stable_fit_calibration.R" else "rmut_stable_fit_validation.R"
  c(
    "rmut_pairing.R" = file.path(common, "rmut_pairing.R"),
    "rmut_stable_fit_helper.R" = file.path(common, stable),
    "data_simulation_final.R" = file.path(Sys.getenv("V33_CODE_ROOT"), "src", "dgm", "primary", case, "data_simulation_final.R"),
    "production_runtime.R" = file.path(common, "production_runtime.R"),
    "rmut_sensitivity_runtime.R" = file.path(common, "rmut_sensitivity_runtime.R")
  )
}

# Runtime support for RMUT sensitivity analyses.
rmut_implementation_version <- "RMUT_COMPLETE_ROW_PAIRING_20261010"

rmut_scenario_code <- function(label) {
  if (grepl("H0", label, fixed = TRUE)) return(1L)
  if (grepl("H1_observed|Realcase", label)) return(3L)
  if (grepl("H1", label)) return(2L)
  # Original IMpower observed scenario is labelled Realcase; reject unknown labels.
  stop("Cannot identify scenario for pairing seed: ", label)
}

rmut_apply_parameter_manifest <- function(case, defaults) {
  path <- Sys.getenv("RMUT_LOCKED_PARAMETER_MANIFEST", unset = "")
  if (!nzchar(path)) {
    return(c(defaults, list(source = "manifest-supplied parameter defaults; no recalibration for this sensitivity run")))
  }
  z <- utils::read.delim(normalizePath(path, mustWork = TRUE), stringsAsFactors = FALSE)
  required <- c("case", "candidate_id", "c", "gamma")
  if (!all(required %in% names(z))) stop("Parameter manifest requires case,candidate_id,c,gamma.")
  if ("method" %in% names(z)) z <- z[z$method == "RMUT", , drop = FALSE]
  z$case[z$case == "IMpower"] <- "IMpower131"
  z$case[z$case == "KEYNOTE598"] <- "KEYNOTE-598"
  required_cases <- c("IMpower131", "KEYNOTE-598", "OAK", "survival_colon")
  if (nrow(z) != 4L || anyDuplicated(z$case) || !setequal(z$case, required_cases)) {
    stop("Parameter manifest must contain exactly one RMUT row for each of the four cases.")
  }
  if (any(!is.finite(z$c)) || any(z$c <= 0 | z$c >= 1) ||
      any(!is.finite(z$gamma)) || any(z$gamma < 0) ||
      any(!is.finite(suppressWarnings(as.numeric(z$candidate_id)))) ||
      any(suppressWarnings(as.numeric(z$candidate_id)) < 1) ||
      any(suppressWarnings(as.numeric(z$candidate_id)) != trunc(suppressWarnings(as.numeric(z$candidate_id))))) stop("Invalid manifest parameters.")
  row <- z[z$case == case, , drop = FALSE]
  list(candidate_id = as.character(row$candidate_id), c = row$c, gamma = row$gamma,
       source = normalizePath(path, mustWork = TRUE), manifest_md5 = unname(tools::md5sum(path)))
}

rmut_preflight <- function(directory, is_rho) {
  files <- c("run_independent_validation_10000.R", "data_simulation_final.R", "rmut_pairing.R",
             "rmut_stable_fit_helper.R", "rmut_sensitivity_runtime.R")
  paths <- c(file.path(directory, files[1L]), unname(rmut_shared_paths(directory)[files[-1L]]))
  if (is_rho) paths <- c(paths, file.path(Sys.getenv("V33_CODE_ROOT"), "src", "common", "rho_parameters.R"))
  if (any(!file.exists(paths))) stop("Missing local dependencies: ", paste(paths[!file.exists(paths)], collapse = ", "))
  for (p in paths) parse(file = p, encoding = "UTF-8")
  if (is_rho) {
    table <- file.path(Sys.getenv("V33_CODE_ROOT"), "config", "rho_parameters.csv")
    if (!file.exists(table)) stop("Missing frozen KM-target parameter table.")
    pe <- new.env(parent = globalenv())
    source(file.path(Sys.getenv("V33_CODE_ROOT"), "src", "common", "rho_parameters.R"), local = pe)
    for (ca in c("IMpower131", "KEYNOTE-598", "OAK", "survival_colon")) {
      for (sc in c("H0", "H1_standard", "H1_observed")) {
        for (arm in c("control", "treatment")) {
          for (lev in c("low", "high")) pe$rho_get_arm(ca, sc, arm, lev)
        }
      }
    }
  }
  invisible(paths)
}

rmut_pp_cache_context <- function(label) {
  files <- c(iv_script_file, iv_dgm_file,
             unname(rmut_shared_paths(iv_version_dir)[c("rmut_pairing.R",
                    "rmut_stable_fit_helper.R", "rmut_sensitivity_runtime.R")]))
  if (exists("iv_parameter_table_path", inherits = TRUE)) {
    files <- c(files, iv_parameter_table_path, file.path(Sys.getenv("V33_CODE_ROOT"), "src", "common", "rho_parameters.R"))
  }
  list(implementation = rmut_implementation_version, analysis = analysis_version,
       case = iv_case, label = label, n_sim = N.sim_realcase, tau = tau_fixed,
       weight = w_R, data_seed = iv_seed_from_label(label), model_seed_base = iv_model_seed_base,
       pairing_seed_policy = "canonical rmut_pair_seed(j,k,scenario_id,case_id)",
       dependency_md5 = unname(tools::md5sum(files)),
       rho = if (exists("iv_sensitivity_rho", inherits = TRUE)) iv_sensitivity_rho else NULL)
}

rmut_cached_stats_valid <- function(stats, expected_context, n, n_looks) {
  if (!is.list(stats) || length(stats) != n ||
      !identical(attr(stats, "rmut_pp_cache_context"), expected_context)) return(FALSE)
  all(vapply(seq_len(n), function(j) {
    x <- stats[[j]]
    is.list(x) && identical(x$replicate_id, as.integer(j)) &&
      is.numeric(x$PP) && length(x$PP) == n_looks &&
      is.numeric(x$PP_final) && length(x$PP_final) == 1L &&
      all((is.na(x$PP) & !is.nan(x$PP)) | (is.finite(x$PP) & x$PP >= 0 & x$PP <= 1)) &&
      ((is.na(x$PP_final) && !is.nan(x$PP_final)) || (is.finite(x$PP_final) && x$PP_final >= 0 && x$PP_final <= 1)) &&
      is.list(x$fit_pair_audit) && length(x$fit_pair_audit) == n_looks + 1L &&
      all(vapply(x$fit_pair_audit, is.list, logical(1)))
  }, logical(1)))
}

rmut_stage_audit <- function(death_fit, composite_fit, pair_seed, tau) {
  da <- attr(death_fit, "rmut_fit_audit")
  ca <- attr(composite_fit, "rmut_fit_audit")
  if (is.null(da) || is.null(ca)) stop("Stable-fit helper did not attach fit audit metadata.")
  list(death = da, composite = ca, pair_seed = as.integer(pair_seed), tau = tau,
       pairing_rule = rmut_implementation_version,
       death_rows = da$posterior_rows, composite_rows = ca$posterior_rows,
       paired_rows = if (isTRUE(da$passed) && isTRUE(ca$passed)) min(da$posterior_rows, ca$posterior_rows) else NA_integer_,
       unmatched_lengths = if (isTRUE(da$passed) && isTRUE(ca$passed)) da$posterior_rows != ca$posterior_rows else NA)
}

rmut_summarize_stats <- function(stats, label, scenario_id, boundaries, interim_counts) {
  n <- length(stats)
  complete <- vapply(stats, function(x) all(is.finite(x$PP)) && is.finite(x$PP_final), logical(1))
  go <- rep(FALSE, n)
  early <- rep(FALSE, n)
  final_no_go <- rep(FALSE, n)
  sample_n <- rep(N.max, n)
  stop_stage <- rep(NA_integer_, n)
  for (j in which(complete)) {
    x <- stats[[j]]
    stages <- which(x$PP <= boundaries)
    if (length(stages)) {
      early[j] <- TRUE
      stop_stage[j] <- stages[1L]
      sample_n[j] <- interim_counts[stages[1L]]
    } else {
      go[j] <- x$PP_final > PP_final_threshold
      final_no_go[j] <- !go[j]
      stop_stage[j] <- length(interim_counts) + 1L
    }
  }
  valid_n <- sum(complete)
  row <- data.frame(
    case = iv_case, scenario = scenario_id, label = label,
    analysis_version = iv_version, parameter_source = iv_parameter_source,
    candidate_id = iv_candidate_id, c = PP_final_threshold, gamma = lambda_prime_bayes,
    total = n, valid = valid_n, fit_failure = n-valid_n, go_count = sum(go),
    go_rate = mean(go), PET = mean(early), final_no_go_rate = mean(final_no_go),
    ESS = mean(sample_n), failure_rule = "no-go; ESS=Nmax; fixed total denominator",
    conditional_valid_go_rate = if (valid_n) mean(go[complete]) else NA_real_,
    conditional_valid_PET = if (valid_n) mean(early[complete]) else NA_real_,
    conditional_valid_ESS = if (valid_n) mean(sample_n[complete]) else NA_real_,
    stringsAsFactors = FALSE)
  replicate <- data.frame(replicate_id = seq_len(n), complete = complete,
                          go = go, early_stop = early, final_no_go = final_no_go,
                          failed_fit = !complete, stop_stage = stop_stage, sample_size = sample_n)
  pp <- do.call(rbind, lapply(stats, function(x) c(x$PP, x$PP_final)))
  colnames(pp) <- c(paste0("PP_", seq_along(interim_counts)), "PP_final")
  replicate <- cbind(replicate, pp)
  utils::write.csv(row, file.path(iv_output_dir, paste0("summary_", scenario_id, ".csv")), row.names = FALSE)
  utils::write.csv(replicate, file.path(iv_output_dir, paste0("replicate_decisions_", scenario_id, ".csv")), row.names = FALSE)
  row
}

rmut_write_completion <- function() {
  if (sink.number(type = "output") != 0L) stop("Close output sink before completion marker.")
  writeLines(c(paste0("version=", iv_version), paste0("implementation=", rmut_implementation_version),
               paste0("case=", iv_case), paste0("run_tag=", iv_run_tag),
               paste0("parameter_source=", iv_parameter_source),
               paste0("completed_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))),
             file.path(iv_output_dir, "RUN_COMPLETE.txt"))
}
