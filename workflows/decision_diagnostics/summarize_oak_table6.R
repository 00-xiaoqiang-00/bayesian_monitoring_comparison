# Reconstruct the OAK H1-observed Table 6 aggregates from saved validation PP.
# Base R only. This script does not simulate patients or fit models.

# R on Windows needs a working UTF-8 locale for the supplied Unicode paths.
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "English_United States.utf8")))
}

args <- commandArgs(trailingOnly = TRUE)
usage <- function() {
  cat(paste(
    "Usage: Rscript summarize_oak_table6.R --cache-dir EXTERNAL_CACHE_DIR --output-dir EXTERNAL_OUTPUT_DIR",
    "       [--parameters current_parameters.csv] [--trial-process trial_process.csv]",
    "       [--schedule four_case_schedule.csv] [--cache-spec oak_table6_cache_spec.csv]",
    "       [--reference archived_OAK_H1observed_table6.csv]",
    "Output directories must be outside this code package. Existing outputs are not overwritten.",
    sep = "\n"), "\n")
}
if (!length(args) || "--help" %in% args || "-h" %in% args) {
  usage()
  quit(status = if (length(args)) 0L else 2L)
}
allowed <- c("--cache-dir", "--output-dir", "--parameters", "--trial-process",
             "--schedule", "--cache-spec", "--reference")
if (length(args) %% 2L != 0L) stop("Options must be supplied as --name value pairs.")
opts <- list()
for (j in seq.int(1L, length(args), by = 2L)) {
  if (!args[j] %in% allowed) stop("Unknown option: ", args[j])
  if (!is.null(opts[[args[j]]])) stop("Duplicate option: ", args[j])
  opts[[args[j]]] <- args[j + 1L]
}
get_opt <- function(name, default = NULL) {
  value <- opts[[name]]
  if (is.null(value)) default else value
}
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- full_args[startsWith(full_args, "--file=")]
if (length(script_arg) != 1L) stop("Run this file through Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
package_root <- dirname(dirname(dirname(script_file)))
parameters_file <- get_opt("--parameters", file.path(package_root, "config/current_parameters.csv"))
trial_file <- get_opt("--trial-process", file.path(package_root, "config/trial_process.csv"))
schedule_file <- get_opt("--schedule", file.path(package_root, "config/analysis_schedule.csv"))
spec_file <- get_opt("--cache-spec", file.path(package_root, "config/oak_table6_cache_spec.csv"))
cache_root <- get_opt("--cache-dir")
output_dir <- get_opt("--output-dir")
if (is.null(cache_root) || is.null(output_dir)) stop("--cache-dir and --output-dir are required.")
cache_root <- normalizePath(cache_root, winslash = "/", mustWork = TRUE)
# Normalize a new output directory through an existing parent to enforce the boundary.
expanded_output <- path.expand(output_dir)
if (!dir.exists(expanded_output)) {
  parent <- dirname(expanded_output)
  if (!dir.exists(parent)) stop("Create the parent directory of --output-dir first.")
  expanded_output <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(expanded_output))
} else {
  expanded_output <- normalizePath(expanded_output, winslash = "/", mustWork = TRUE)
}
inside <- function(path, root) {
  path <- tolower(gsub("\\\\", "/", path))
  root <- tolower(gsub("\\\\", "/", root))
  identical(path, root) || startsWith(path, paste0(root, "/"))
}
submission_root <- dirname(package_root)
protected_root <- if (dir.exists(file.path(submission_root, "04_Data_and_Results"))) submission_root else package_root
if (inside(expanded_output, protected_root)) stop("--output-dir must be outside the entire submission (or standalone code package).")
output_names <- c("OAK_H1observed_table6.csv", "OAK_H1observed_PP_quantiles.csv", "input_provenance.csv")
if (any(file.exists(file.path(expanded_output, output_names)))) stop("An output file already exists; choose a new output directory.")

parameters <- read.csv(parameters_file, stringsAsFactors = FALSE, check.names = FALSE)
trial <- read.csv(trial_file, stringsAsFactors = FALSE, check.names = FALSE)
schedule <- read.csv(schedule_file, stringsAsFactors = FALSE, check.names = FALSE)
spec <- read.csv(spec_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(as.character(spec$method), c("Win Ratio", "Cox", "RMST", "RMUT"))) stop("Cache specification must list Win Ratio, Cox, RMST, RMUT in that order.")
if (!all(c("relative_cache", "format") %in% names(spec))) stop("Invalid cache specification.")
oak_trial <- trial[trial$case == "OAK", , drop = FALSE]
if (nrow(oak_trial) != 1L) stop("Expected exactly one OAK trial row.")
oak_schedule <- schedule[schedule$case == "OAK" & grepl("主分析|最终分析", schedule$analysis), , drop = FALSE]
if (nrow(oak_schedule) != 4L) stop("Expected three OAK main interim looks and one final look.")
cohorts <- as.numeric(oak_schedule$trigger_n)
if (any(!is.finite(cohorts)) || !all(diff(cohorts) > 0) || cohorts[4] != oak_trial$Nmax) stop("Invalid OAK look cohorts.")
fractions <- cohorts[1:3] / oak_trial$Nmax
if (!isTRUE(all.equal(fractions, as.numeric(strsplit(oak_trial$interim_information_fractions, ";", fixed = TRUE)[[1]]), tolerance = 1e-12))) stop("Trial and schedule information fractions disagree.")

rows <- list()
quantile_rows <- list()
input_paths <- c(script = script_file, parameters = parameters_file, trial_process = trial_file,
                 schedule = schedule_file, cache_spec = spec_file)
for (i in seq_len(nrow(spec))) {
  method <- spec$method[i]
  parameter <- parameters[parameters$case == "OAK" & parameters$method == method, , drop = FALSE]
  if (nrow(parameter) != 1L) stop("Expected one OAK threshold row for ", method)
  eta <- parameter$c[[1]]
  gamma <- parameter$gamma[[1]]
  if (!is.finite(eta) || !is.finite(gamma)) stop("Non-finite threshold.")
  cache_file <- file.path(cache_root, spec$relative_cache[i])
  input_paths[paste0("cache_", gsub(" ", "_", method))] <- cache_file
  cache <- readRDS(cache_file)
  if (spec$format[i] == "winratio_matrix") {
    if (!is.matrix(cache$PP) || ncol(cache$PP) != 4L) stop("Expected four-look WR PP matrix.")
    all_pp <- cache$PP
    fit_ok <- rep(TRUE, nrow(all_pp))
    if (!is.null(cache$look_cohort) && !isTRUE(all.equal(as.numeric(cache$look_cohort), cohorts, tolerance = 1e-12))) stop("WR cohort schedule differs from locked schedule.")
    if (!is.null(cache$N.max) && cache$N.max != oak_trial$Nmax) stop("WR N.max differs from locked trial.")
  } else if (spec$format[i] == "stats_list") {
    if (!is.list(cache) || !length(cache)) stop("Expected a nonempty stats list.")
    all_pp <- t(vapply(cache, function(entry) {
      if (length(entry$PP) != 3L || length(entry$PP_final) != 1L) stop("Expected three interim PP and one PP_final per replicate.")
      as.numeric(c(entry$PP, entry$PP_final))
    }, numeric(4L)))
    fit_ok <- vapply(cache, function(entry) {
      if (!is.null(entry$ok)) return(isTRUE(entry$ok))
      if (!is.null(entry$fit_ok)) return(isTRUE(entry$fit_ok))
      TRUE
    }, logical(1))
  } else stop("Unknown cache format: ", spec$format[i])
  valid <- fit_ok & rowSums(is.finite(all_pp)) == 4L
  if (!any(valid)) stop("No valid replicates for ", method)
  pp <- all_pp[valid, , drop = FALSE]
  if (any(pp < 0 | pp > 1)) stop("PP outside [0,1].")
  boundaries <- eta * fractions ^ gamma
  stop1 <- pp[, 1] <= boundaries[1]
  stop2 <- !stop1 & pp[, 2] <= boundaries[2]
  stop3 <- !stop1 & !stop2 & pp[, 3] <= boundaries[3]
  final <- !(stop1 | stop2 | stop3)
  # Exact rule of the locked strategy and manuscript: strict > at final.
  final_success <- pp[, 4] > eta
  go <- final & final_success
  n <- nrow(pp)
  stop_counts <- c(sum(stop1), sum(stop2), sum(stop3))
  if (sum(stop_counts) + sum(final) != n) stop("Path decomposition failed.")
  rows[[i]] <- data.frame(
    case = "OAK", scenario = "H1_observed", method = method,
    n_stored = nrow(all_pp), n_valid = n, n_invalid = sum(!valid),
    eta = eta, gamma = gamma, b1 = boundaries[1], b2 = boundaries[2], b3 = boundaries[3],
    pp1_median = median(pp[, 1]), pp2_median = median(pp[, 2]),
    pp3_median = median(pp[, 3]), pp_final_median = median(pp[, 4]),
    stop1 = stop_counts[1], stop2 = stop_counts[2], stop3 = stop_counts[3],
    final = sum(final), go = sum(go), final_no_go = sum(final & !final_success),
    PET = sum(stop_counts) / n, go_probability = mean(go),
    ESS = mean(ifelse(stop1, cohorts[1], ifelse(stop2, cohorts[2], ifelse(stop3, cohorts[3], cohorts[4])))),
    diagnostic_final_gt_count = sum(final_success), diagnostic_final_gt_probability = mean(final_success),
    final_equal_eta_count = sum(pp[, 4] == eta),
    reached_final_equal_eta_count = sum(final & pp[, 4] == eta),
    stopped1_but_final_gt_probability = if (any(stop1)) mean(final_success[stop1]) else NA_real_,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  quantile_rows[[i]] <- do.call(rbind, lapply(seq_len(4L), function(k) data.frame(
    case = "OAK", scenario = "H1_observed", method = method, look = k,
    quantile = c(0.05, 0.25, 0.50, 0.75, 0.95),
    PP = as.numeric(quantile(pp[, k], c(0.05, 0.25, 0.50, 0.75, 0.95), type = 7)),
    stringsAsFactors = FALSE
  )))
}
summary <- do.call(rbind, rows)
if (!is.null(get_opt("--reference"))) {
  reference <- read.csv(get_opt("--reference"), stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(reference), names(summary)) || !identical(as.character(reference$method), as.character(summary$method))) stop("Reference columns or method order differ.")
  for (column in names(summary)) {
    if (is.numeric(summary[[column]])) {
      if (!isTRUE(all.equal(reference[[column]], summary[[column]], tolerance = 1e-12, check.attributes = FALSE))) stop("Reference differs in column ", column)
    } else if (!identical(as.character(reference[[column]]), as.character(summary[[column]]))) stop("Reference differs in column ", column)
  }
  cat("Reference check passed.\n")
}
if (!dir.exists(expanded_output)) dir.create(expanded_output, recursive = FALSE)
write.csv(summary, file.path(expanded_output, output_names[1]), row.names = FALSE, na = "")
write.csv(do.call(rbind, quantile_rows), file.path(expanded_output, output_names[2]), row.names = FALSE, na = "")
provenance <- data.frame(
  role = names(input_paths),
  path = vapply(input_paths, normalizePath, character(1), winslash = "/", mustWork = TRUE),
  md5 = as.character(unname(tools::md5sum(input_paths))), stringsAsFactors = FALSE
)
write.csv(provenance, file.path(expanded_output, output_names[3]), row.names = FALSE)
print(summary[, c("method", "n_valid", "stop1", "stop2", "stop3", "final", "go", "PET", "diagnostic_final_gt_probability", "final_equal_eta_count", "reached_final_equal_eta_count")], row.names = FALSE, digits = 8)
cat("Results written to:", expanded_output, "\n")
