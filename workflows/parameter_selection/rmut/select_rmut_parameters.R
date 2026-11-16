#!/usr/bin/env Rscript
# Verify a saved formal selection, or explicitly rescreen complete calibration PP caches.
# No patient-data generator or model-fitting code is loaded.

options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))

selection_details <- function(screened) {
  eligible <- screened$eligible
  design <- screened$design
  h0 <- screened$h0
  strict <- eligible[which.min(eligible$H0_ESS), , drop = FALSE]
  n_ref <- n_used_vector(h0, strict$c, strict$gamma, design$interim, design$n_max)
  se <- vapply(seq_len(nrow(eligible)), function(i) {
    n_i <- n_used_vector(h0, eligible$c[i], eligible$gamma[i], design$interim, design$n_max)
    stats::sd(n_i - n_ref) / sqrt(h0$attempted)
  }, numeric(1L))
  eligible$strict_min_candidate_id <- strict$candidate_id
  eligible$delta_H0_ESS <- eligible$H0_ESS - strict$H0_ESS
  eligible$paired_SE_delta <- se
  eligible$near_optimal_2SE <- eligible$H0_ESS <= strict$H0_ESS + 2 * se
  eligible$grid_boundary <- eligible$c <= min(pp_grid) + 1e-9 |
    eligible$c >= max(pp_grid) - 1e-9 |
    eligible$gamma <= 0.01 + 1e-9 |
    eligible$gamma >= max(lambda_grid) - 1e-9
  eligible$legacy_boundary_flag <- eligible$c <= 0.855 | eligible$gamma <= 0.015
  eligible
}

read_cache_manifest <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  required <- c("case", "scenario", "cache_path", "expected_attempted", "sha256", "md5")
  if (!all(required %in% names(manifest))) stop("Incomplete cache manifest")
  keys <- paste(manifest$case, manifest$scenario)
  expected <- as.vector(outer(names(designs), c("H0", "H1-standard"), paste))
  if (anyDuplicated(keys) || !setequal(keys, expected) || nrow(manifest) != 8L) {
    stop("Cache manifest must contain exactly H0 and H1-standard for all four cases")
  }
  if (any(manifest$expected_attempted != 10000L)) stop("Each formal cache must specify 10,000 attempts")
  manifest$absolute_path <- vapply(manifest$cache_path, function(p) {
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", p)) p <- file.path(dirname(path), p)
    normalizePath(p, winslash = "/", mustWork = TRUE)
  }, character(1L))
  actual <- unname(tools::md5sum(manifest$absolute_path))
  if (anyNA(actual) || any(tolower(actual) != tolower(manifest$md5))) {
    stop("Cache bytes do not match the manifest MD5; SHA-256 identities are also recorded in the manifest")
  }
  manifest
}

verify_saved_case <- function(case, selection_dir) {
  eligible_path <- file.path(selection_dir, paste0("eligible_", case, ".tsv"))
  selected_path <- file.path(selection_dir, paste0("selected_", case, ".tsv"))
  eligible <- utils::read.delim(eligible_path, check.names = FALSE)
  saved <- utils::read.delim(selected_path, check.names = FALSE)
  required <- c("case", "candidate_id", "c", "gamma", "H0_T1E", "H0_PET", "H0_ESS", "H1_Power", "H1_PET", "H1_ESS")
  if (!all(required %in% names(eligible)) || !nrow(eligible) ||
      anyDuplicated(eligible$candidate_id) || any(eligible$case != case) ||
      nrow(saved) != 1L || !identical(as.character(saved$case), case)) stop("Invalid saved selection table: ", case)
  if (any(!is.finite(as.matrix(eligible[, setdiff(required, "case"), drop = FALSE])))) {
    stop("Non-finite saved selection values: ", case)
  }
  i <- match(eligible$c, pp_grid)
  j <- match(eligible$gamma, lambda_grid)
  if (anyNA(i) || anyNA(j) || any(eligible$candidate_id != i + (j - 1L) * length(pp_grid))) {
    stop("Candidate IDs do not match the complete grid: ", case)
  }
  if (any(eligible$H0_T1E != 0.05 | eligible$H1_Power < 0.80 |
          eligible$c < 0.85 | eligible$gamma < 0.01)) stop("Saved eligible rows fail a selection gate: ", case)
  h0 <- prepare_cache(cache_path(case, "H0"))
  h1 <- prepare_cache(cache_path(case, "H1"))
  d <- designs[[case]]
  max_error <- 0
  for (k in seq_len(nrow(eligible))) {
    m0 <- scenario_metrics(h0, eligible$c[k], eligible$gamma[k], d$interim, d$n_max)
    m1 <- scenario_metrics(h1, eligible$c[k], eligible$gamma[k], d$interim, d$n_max)
    if (m0[["success_n"]] != 500L) stop("Saved candidate no longer has exactly 500 H0 go decisions: ", case)
    stored <- as.numeric(eligible[k, c("H0_T1E", "H0_PET", "H0_ESS", "H1_Power", "H1_PET", "H1_ESS")])
    recomputed <- unname(c(m0[c("rate", "pet", "ess")], m1[c("rate", "pet", "ess")]))
    max_error <- max(max_error, abs(stored - recomputed))
  }
  if (max_error > 1e-8) stop("Saved eligible metrics differ from PP-cache classification: ", case)
  screened <- list(eligible = eligible, h0 = h0, design = d)
  chosen <- select_case(screened)
  for (field in intersect(names(saved), names(chosen))) {
    if (is.numeric(chosen[[field]])) {
      if (any(abs(saved[[field]] - chosen[[field]]) > 1e-8)) stop("Selected field differs: ", case, " / ", field)
    } else if (!identical(as.character(saved[[field]]), as.character(chosen[[field]]))) {
      stop("Selected field differs: ", case, " / ", field)
    }
  }
  details <- selection_details(screened)
  row <- data.frame(
    case = case, eligible_n = nrow(eligible),
    strict_min_count = sum(eligible$H0_ESS == min(eligible$H0_ESS)),
    strict_min_candidate_id = details$strict_min_candidate_id[1L],
    selected_candidate_id = chosen$candidate_id,
    c = chosen$c, gamma = chosen$gamma,
    near_n = sum(details$near_optimal_2SE),
    grid_boundary = chosen$boundary,
    legacy_boundary_flag = chosen$c <= 0.855 | chosen$gamma <= 0.015,
    H0_valid = sum(h0$valid), H1_valid = sum(h1$valid),
    attempted_per_scenario = h0$attempted,
    max_metric_error = max_error,
    selected_matches_saved = TRUE,
    stringsAsFactors = FALSE
  )
  list(summary = row, candidates = details)
}

new_output_directory <- function(path, code_root) {
  source(file.path(code_root,"src","common","output_paths.R"),local=TRUE)
  path <- external_output_path(path,code_root)
  if (!nzchar(path) || file.exists(path) || dir.exists(path)) stop("Output must name a new, nonexistent directory")
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_selection_metadata <- function(out_dir, mode, manifest) {
  lines <- c(
    "selection_rule=FORMAL_RMUT_ROWPAIR_20261010",
    paste0("operation=", mode),
    "input=complete_H0_and_H1_standard_calibration_PP_paths",
    "patient_data_generation=false", "model_fitting=false",
    "attempted_per_scenario=10000", "failed_fit_rule=no_go_and_Nmax",
    "required_valid_rate_each_scenario=0.95",
    "grid_c=0.80_to_0.99_by_0.0005", "grid_gamma=0_to_5_by_0.005",
    "candidate_id=c_index_plus_gamma_index_minus_one_times_381",
    "H0_go_count=500", "H1_power_min=0.80", "c_min=0.85", "gamma_min=0.01",
    "strict_reference=first_minimum_H0_ESS_in_original_grid_order",
    "near_rule=H0_ESS_delta_le_2_paired_MCSE",
    "near_ranking=H1_Power_desc,c_desc,H1_PET_asc,grid_boundary_asc,candidate_id_asc",
    "grid_boundary=c_le_0.80_or_c_ge_0.99_or_gamma_le_0.01_or_gamma_ge_5_with_1e-9_tolerance",
    "legacy_boundary_flag=c_le_0.855_or_gamma_le_0.015_diagnostic_only",
    "H1_observed_used=false",
    if (mode == "verify_saved_selection") "complete_grid_rescreen_executed=false" else "complete_grid_rescreen_executed=true"
  )
  writeLines(lines, file.path(out_dir, "SELECTION_METADATA.txt"), useBytes = TRUE)
  utils::write.csv(manifest[, setdiff(names(manifest), "absolute_path"), drop = FALSE],
                   file.path(out_dir, "INPUT_CACHE_IDENTITIES.csv"), row.names = FALSE)
}

rmut_selection_usage <- function() {
  cat(c(
    "RMUT parameter selection from fixed calibration PP inputs.",
    "",
    "Usage:",
    "  Rscript workflows/parameter_selection/rmut/select_rmut_parameters.R --help",
    "  Rscript workflows/parameter_selection/rmut/select_rmut_parameters.R rescreen NEW_OUTPUT_DIR [CASE ...]",
    "  Rscript workflows/parameter_selection/rmut/select_rmut_parameters.R verify --reference-dir DIR [--output-dir NEW_OUTPUT_DIR]",
    "",
    "Cases: IMpower, KEYNOTE598, OAK, survival_colon. Omit CASE to rescreen all four.",
    "verify requires eligible_<case>.tsv and selected_<case>.tsv for all four cases.",
    "Reference paths resolve from the current working directory; no submission layout is required.",
    "Output directories must be new and outside this repository.",
    "No arguments or --help shows this help without reading the PP inputs.",
    ""
  ), sep = "\n")
  invisible(TRUE)
}

parse_verify_arguments <- function(args) {
  values <- list(reference_dir = NULL, output_dir = NULL)
  positional <- character()
  flags <- c("--reference-dir" = "reference_dir", "--output-dir" = "output_dir")
  set_value <- function(name, value) {
    if (!is.null(values[[name]])) stop("Duplicate option: ", gsub("_", "-", name), call. = FALSE)
    if (!nzchar(value) || startsWith(value, "--")) stop("Missing path for --", gsub("_", "-", name), call. = FALSE)
    values[[name]] <<- value
  }
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% names(flags)) {
      if (i == length(args)) stop("Missing path for ", arg, call. = FALSE)
      set_value(unname(flags[[arg]]), args[[i + 1L]])
      i <- i + 2L
      next
    }
    equals_flag <- names(flags)[vapply(names(flags), function(flag) startsWith(arg, paste0(flag, "=")), logical(1L))]
    if (length(equals_flag)) {
      flag <- equals_flag[[1L]]
      set_value(unname(flags[[flag]]), substring(arg, nchar(flag) + 2L))
    } else if (startsWith(arg, "-")) {
      stop("Unknown verify option: ", arg, ". Use --help for usage.", call. = FALSE)
    } else {
      positional <- c(positional, arg)
    }
    i <- i + 1L
  }
  if (length(positional) > 1L) stop("verify accepts at most one positional output directory", call. = FALSE)
  if (length(positional)) set_value("output_dir", positional[[1L]])
  if (is.null(values$reference_dir)) stop("verify requires --reference-dir DIR containing the eight eligible/selected tables. Use rescreen to generate them.", call. = FALSE)
  values
}

read_reference_directory <- function(path, cases) {
  path <- path.expand(path)
  if (!dir.exists(path)) stop("Reference directory does not exist: ", path, call. = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  required <- as.vector(rbind(paste0("eligible_", cases, ".tsv"), paste0("selected_", cases, ".tsv")))
  files <- file.path(path, required)
  readable <- file.exists(files) & !dir.exists(files) & file.access(files, 4L) == 0L
  if (any(!readable)) stop("Reference directory must contain eight readable tables. Missing or unreadable: ", paste(required[!readable], collapse = ", "), call. = FALSE)
  path
}

rmut_selection_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || any(args %in% c("--help", "-h", "help"))) {
    rmut_selection_usage()
    return(invisible(TRUE))
  }
  mode <- args[[1L]]
  if (!mode %in% c("verify", "rescreen")) stop("Unknown mode: ", mode, ". Use --help for usage.", call. = FALSE)
  verify_args <- if (mode == "verify") parse_verify_arguments(args[-1L]) else NULL
  if (mode == "rescreen" && length(args) < 2L) stop("rescreen requires a new output directory. Use --help for usage.", call. = FALSE)
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Run this entry with Rscript")
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
  source(file.path(script_dir, "..", "..", "..", "src", "common", "parameter_selection", "rmut_selection_functions.R"), local = .GlobalEnv, encoding = "UTF-8")
  if (mode == "verify") selection_dir <- read_reference_directory(verify_args$reference_dir, names(designs))
  default_manifest <- file.path(script_dir, "..", "..", "..", "inputs", "calibration_pp", "input_manifest.csv")
  manifest_path <- Sys.getenv("RMUT_FIXED_INPUT_MANIFEST", unset = default_manifest)
  manifest <- read_cache_manifest(manifest_path)
  cache_path <<- function(case, kind) {
    scenario <- if (kind == "H0") "H0" else "H1-standard"
    hit <- manifest$absolute_path[manifest$case == case & manifest$scenario == scenario]
    if (length(hit) != 1L) stop("Cache lookup failed: ", case, " / ", scenario)
    hit
  }
  if (mode == "verify") {
    checks <- lapply(names(designs), verify_saved_case, selection_dir = selection_dir)
    summary <- do.call(rbind, lapply(checks, `[[`, "summary"))
    print(summary, row.names = FALSE)
    if (!is.null(verify_args$output_dir)) {
      out_dir <- new_output_directory(verify_args$output_dir, file.path(script_dir,"..","..",".."))
      utils::write.csv(summary, file.path(out_dir, "saved_selection_verification.csv"), row.names = FALSE)
      utils::write.table(do.call(rbind, lapply(checks, `[[`, "candidates")),
                         file.path(out_dir, "eligible_with_paired_mcse.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      write_selection_metadata(out_dir, "verify_saved_selection", manifest)
    }
    message("SAVED_SELECTION_VERIFIED; complete-grid rescreen was not requested")
  } else {
    if (length(args) < 2L) stop("rescreen requires a new output directory")
    requested <- if (length(args) >= 3L) args[3L:length(args)] else names(designs)
    if (anyDuplicated(requested) || any(!requested %in% names(designs))) stop("Name each requested case once")
    out_dir <- new_output_directory(args[[2L]], file.path(script_dir,"..","..",".."))
    chosen_rows <- lapply(requested, function(case) {
      screened <- screen_case(case)
      chosen <- select_case(screened)
      utils::write.table(screened$eligible, file.path(out_dir, paste0("eligible_", case, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      utils::write.table(selection_details(screened), file.path(out_dir, paste0("eligible_with_paired_mcse_", case, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      utils::write.table(chosen, file.path(out_dir, paste0("selected_", case, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      chosen
    })
    selected <- do.call(rbind, chosen_rows)
    manifest_name <- if (setequal(requested, names(designs))) "RMUT_LOCKED_PARAMETERS.tsv" else "SELECTED_PARAMETER_SUBSET.tsv"
    utils::write.table(selected[, c("case", "candidate_id", "c", "gamma", "status")], file.path(out_dir, manifest_name), sep = "\t", quote = FALSE, row.names = FALSE)
    write_selection_metadata(out_dir, "complete_grid_rescreen", manifest)
    print(selected, row.names = FALSE)
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) rmut_selection_main()
