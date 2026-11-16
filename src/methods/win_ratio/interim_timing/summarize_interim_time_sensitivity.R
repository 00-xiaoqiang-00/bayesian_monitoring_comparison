#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = "") {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Pass at most one ", prefix, " argument.")
  if (!length(hit)) return(default)
  sub(paste0("^", prefix), "", hit[1L])
}
script_hit <- commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))]
script_dir <- if (length(script_hit) == 1L) normalizePath(dirname(sub("^--file=", "", script_hit[1L])), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
mode_hit <- args[grepl("^--mode=", args)]
run_mode <- if (length(mode_hit) == 1L) {
  sub("^--mode=", "", mode_hit[1L])
} else {
  "formal_10000"
}

sens_root <- normalizePath(arg_value("--sensitivity-root=", Sys.getenv("SENSITIVITY_ROOT", unset = script_dir)), winslash = "/", mustWork = TRUE)
source(file.path(sens_root, "R", "interim_sensitivity_core.R"))

output_base <- arg_value("--output-root=", Sys.getenv("WINRATIO_S8_OUTPUT_ROOT", unset = ""))
if (!nzchar(output_base)) stop("Pass --output-root=/external/directory or set WINRATIO_S8_OUTPUT_ROOT.")
cases <- c("IMpower131", "KEYNOTE598", "OAK", "survival_colon")
oc_list <- list()
delta_list <- list()
for (case_key in cases) {
  oc_path <- file.path(output_base, run_mode, case_key, "results", "oc_latest.csv")
  delta_path <- file.path(output_base, run_mode, case_key, "results",
                          "paired_delta_latest.csv")
  if (file.exists(oc_path)) {
    oc_list[[case_key]] <- utils::read.csv(oc_path, stringsAsFactors = FALSE)
  }
  if (file.exists(delta_path)) {
    delta_list[[case_key]] <- utils::read.csv(delta_path, stringsAsFactors = FALSE)
  }
}

out_dir <- file.path(output_base, run_mode, "_summary")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (length(oc_list)) {
  oc <- do.call(rbind, oc_list)
  utils::write.csv(oc, file.path(out_dir, "oc_all_cases.csv"), row.names = FALSE)
  print(oc)
}
if (length(delta_list)) {
  delta <- do.call(rbind, delta_list)
  delta$primary_analysis <- delta$scenario != "H1_observed"
  delta$note <- ifelse(
    delta$scenario == "H1_observed",
    "descriptive; do not retune c_F/gamma_E",
    "primary locked-design sensitivity"
  )
  utils::write.csv(delta, file.path(out_dir, "paired_delta_all_cases.csv"),
                   row.names = FALSE)
  print(delta)
}

cat("Summary written to ", out_dir, "\n", sep = "")
