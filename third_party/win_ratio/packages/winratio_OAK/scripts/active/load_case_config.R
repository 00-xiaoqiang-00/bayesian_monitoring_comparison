# Load a case config and prepare isolated output directories.
#
# Usage:
#   source("scripts/active/load_case_config.R")
#   cfg <- load_case_config("OAK_A_ITT")
#   # or: cfg <- load_case_config(commandArgs(TRUE))

load_case_config <- function(case_id = NULL,
                             repo_root = getwd(),
                             args = commandArgs(TRUE)) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  if (is.null(case_id) || !nzchar(case_id)) {
    # support: --case=OAK_A_ITT  or  case=OAK_A_ITT  or  first bare arg
    m <- regexec("(?:--)?case=([^[:space:]]+)", paste(args, collapse = " "))
    hit <- regmatches(paste(args, collapse = " "), m)[[1]]
    if (length(hit) >= 2) {
      case_id <- hit[2]
    } else if (length(args) >= 1 && !grepl("=", args[1])) {
      case_id <- args[1]
    } else {
      stop("No case_id. Pass --case=OAK_A_ITT or load_case_config(\"OAK_A_ITT\").")
    }
  }

  config_path <- file.path(repo_root, "cases", case_id, "config.R")
  if (!file.exists(config_path)) {
    stop("Missing case config: ", config_path)
  }

  env <- new.env(parent = baseenv())
  sys.source(config_path, envir = env)

  if (is.null(env$case_id) || !identical(as.character(env$case_id), as.character(case_id))) {
    warning("config case_id=\"", env$case_id, "\" differs from requested \"", case_id, "\". Using folder name.")
    env$case_id <- case_id
  }

  output_root <- Sys.getenv("WINRATIO_OUTPUT_ROOT", unset = file.path(repo_root, "output"))
  case_root <- file.path(output_root, "cases", case_id)
  env$repo_root <- repo_root
  env$case_root <- case_root
  env$config_path <- config_path
  env$out_dir <- case_root
  env$cache_dir <- file.path(case_root, "cache")
  env$latest_dir <- file.path(case_root, "latest")
  env$history_dir <- file.path(case_root, "history")
  env$validation_dir <- file.path(case_root, "validation")
  env$logs_dir <- file.path(case_root, "logs")
  env$runs_dir <- file.path(case_root, "runs")

  for (d in c(env$case_root, env$cache_dir, env$latest_dir, env$history_dir,
              env$validation_dir, env$logs_dir, env$runs_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # Soft validation: warn on incomplete placeholders
  na_fields <- c()
  for (nm in c("N.max", "recruit.int", "cens_upper")) {
    if (is.null(env[[nm]]) || (length(env[[nm]]) == 1 && is.na(env[[nm]]))) {
      na_fields <- c(na_fields, nm)
    }
  }
  if (any(is.na(env$mu_h0_trt)) || any(is.na(env$mu_h1_trt)) ||
      any(is.na(env$Sigma.ctrl)) || any(is.na(env$Sigma_h1_trt))) {
    na_fields <- c(na_fields, "mu/Sigma")
  }
  if (identical(env$generator_type, "oak_lognormal_mixture_v1") &&
      is.null(env$oak_generator_parameters)) {
    na_fields <- c(na_fields, "oak_generator_parameters")
  }
  if (is.null(env$final_analysis_time) ||
      !is.finite(env$final_analysis_time)) {
    na_fields <- c(na_fields, "final_analysis_time")
  }
  if (length(na_fields)) {
    message("Case ", case_id, ": incomplete config fields: ",
            paste(unique(na_fields), collapse = ", "))
  }

  cat("Loaded case:", case_id, "\n")
  cat("  config:", config_path, "\n")
  cat("  output:", case_root, "\n")
  env
}

# Convenience: list available cases
list_cases <- function(repo_root = getwd()) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  cases_dir <- file.path(repo_root, "cases")
  if (!dir.exists(cases_dir)) return(character(0))
  names <- list.dirs(cases_dir, full.names = FALSE, recursive = FALSE)
  names[file.exists(file.path(cases_dir, names, "config.R"))]
}


