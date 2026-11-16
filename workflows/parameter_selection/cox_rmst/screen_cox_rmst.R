#!/usr/bin/env Rscript
# Formal exact-count screening.
options(stringsAsFactors = FALSE)
entry_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(entry_arg) != 1L) stop("Run with Rscript.")
entry_dir <- dirname(normalizePath(sub("^--file=", "", entry_arg), mustWork = TRUE))
source(file.path(entry_dir, "..", "..", "..", "src", "common", "parameter_selection", "selection_io.R"), local = TRUE, encoding = "UTF-8")

main <- function() {
  a <- selector_args("Rscript screen_cox_rmst.R --cache-manifest CACHE_MAP.tsv --output-dir NEW_EXTERNAL_DIRECTORY [--method Cox|RMST] [--case CASE]")
  mapping <- selector_manifest(a[["cache-manifest"]], a[["method"]], a[["case"]])
  project_info <- mapping$project_info
  cache_table <- mapping$cache_table
  case_interims <- mapping$case_interims
  source(file.path(entry_dir, "..", "..", "..", "src", "common", "parameter_selection", "selection_functions.R"), local = TRUE, encoding = "UTF-8")
  output_dir <- selector_new_output(a[["output-dir"]], file.path(entry_dir, "..", "..", ".."))
  projects <- unique(project_info$project)
  # Sequential execution is portable to Windows and does not change the grid.
  sets <- lapply(projects, function(nm) {
    message("Exact T1E screening ", nm)
    screen_one(nm)
  })
  names(sets) <- projects
  eligible <- do.call(rbind, sets)
  eligible <- eligible[order(eligible$case, eligible$method, eligible$PP_final, eligible$lambda_prime), , drop = FALSE]
  utils::write.table(eligible, file.path(output_dir, "eligible_parameter_groups.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  metadata <- c(
    "algorithm=NONWIN_EXACT_EMPIRICAL_T1E5_POWER80_V1_20260806;scope=Cox_RMST",
    "H0_target_go_count=500;H0_denominator=10000_attempted",
    "H1_power_min=0.80_on_all_attempted;valid_fit_rate_min=0.95_each_scenario",
    "failure_rule=no_go_and_Nmax_for_ESS;H1_observed_used=false",
    "Cox_grid=c_0.80_to_0.99_by_0.001;gamma_0_to_5_by_0.005",
    "RMST_grid=c_0.80_to_0.99_by_0.001;gamma_0_to_15_by_0.01",
    "status=ELIGIBLE_SET_NOT_RANKED_NOT_LOCKED",
    paste0("input_manifest=", mapping$manifest_path),
    paste0("input_manifest_md5=", unname(tools::md5sum(mapping$manifest_path))),
    paste0("eligible_rows=", nrow(eligible)))
  writeLines(metadata, file.path(output_dir, "SCREENING_METADATA.txt"))
  cat("SCREENING_COMPLETE projects=", length(projects), " eligible_rows=", nrow(eligible), "\n", sep = "")
}
if (sys.nframe() == 0L) main()
