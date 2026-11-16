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
v33_run("src/methods/win_ratio/validation.R","Win Ratio","survival_colon","independent validation",c("--case=survival_colon"))
