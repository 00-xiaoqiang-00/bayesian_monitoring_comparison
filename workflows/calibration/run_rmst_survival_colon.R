# Run this analysis with --output-dir outside the submission, or --preflight-only.
.workflow_hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(.workflow_hit) != 1L) stop("Run with Rscript.")
.workflow_root <- dirname(normalizePath(sub("^--file=", "", .workflow_hit), winslash = "/", mustWork = TRUE))
while (!file.exists(file.path(.workflow_root, "config/current_parameters.csv"))) {
  .workflow_parent <- dirname(.workflow_root)
  if (identical(.workflow_parent, .workflow_root)) stop("Cannot locate the code package.")
  .workflow_root <- .workflow_parent
}
Sys.setenv(V33_CODE_ROOT = .workflow_root)
source(file.path(.workflow_root, "src/common/run_workflow.R"))
v33_run("src/methods/rmst/calibration/survival_colon/calibration/build_calibration_cache_10000.R","RMST","survival_colon","calibration",character())
