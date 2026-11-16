invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
args <- commandArgs(trailingOnly=TRUE)
value <- function(key, default=NULL) {
  i <- match(key,args)
  if(is.na(i)) { if(is.null(default)) stop("Missing ",key); return(default) }
  if(i==length(args)) stop("Missing value for ",key)
  args[i+1L]
}
script_arg <- commandArgs(trailingOnly=FALSE)[grepl("^--file=",commandArgs(trailingOnly=FALSE))]
if(length(script_arg)!=1L) stop("Run this entry with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=","",script_arg[1L]),winslash="/",mustWork=TRUE))
code_root <- normalizePath(value("--code-root",file.path(script_dir,"..","..","..")),winslash="/",mustWork=TRUE)
source(file.path(code_root,"src","common","output_paths.R"),local=TRUE)
output_dir <- external_output_path(value("--output-dir"),code_root)

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
input_file <- file.path(code_root,"inputs","clinical_data","survival_colon.csv")
source(file.path(code_root,"src","dgm","colon_input.R"),local=TRUE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(MASS))

set.seed(20260721)

n_per_arm <- 300000L
target_hr_recurrence <- 0.59
target_hr_death <- 0.67

mu_ctrl <- log(c(38.1268, 64.5907))
sigma_common <- matrix(
  c(3.29767175, 2.19701400,
    2.19701400, 1.77436755),
  nrow = 2, byrow = TRUE
)

z_ctrl <- MASS::mvrnorm(n_per_arm, mu = c(0, 0), Sigma = sigma_common)
z_trt <- MASS::mvrnorm(n_per_arm, mu = c(0, 0), Sigma = sigma_common)

cox_log_hr <- function(time_ctrl, status_ctrl, time_trt, status_trt) {
  time <- c(time_ctrl, time_trt)
  status <- c(status_ctrl, status_trt)
  group <- c(rep.int(0L, length(time_ctrl)), rep.int(1L, length(time_trt)))
  ord <- order(time)
  status <- status[ord]
  group <- group[ord]

  risk_0 <- rev(cumsum(rev(group == 0L)))
  risk_1 <- rev(cumsum(rev(group == 1L)))
  event <- status == 1L

  score <- function(beta) {
    eb <- exp(beta)
    sum(group[event] - eb * risk_1[event] /
          (risk_0[event] + eb * risk_1[event]))
  }

  uniroot(score, interval = c(-4, 4), tol = 1e-10)$root
}

log_ctrl <- sweep(z_ctrl, 2, mu_ctrl, "+")
death_ctrl <- exp(log_ctrl[, 2])
recurrence_ctrl <- exp(log_ctrl[, 1])
recurrence_status_ctrl <- as.integer(recurrence_ctrl <= death_ctrl)
recurrence_followup_ctrl <- pmin(recurrence_ctrl, death_ctrl)

death_log_hr_for_shift <- function(delta_d) {
  death_trt <- exp(mu_ctrl[2] + delta_d + z_trt[, 2])
  cox_log_hr(
    death_ctrl, rep.int(1L, n_per_arm),
    death_trt, rep.int(1L, n_per_arm)
  )
}

delta_d <- uniroot(
  function(x) death_log_hr_for_shift(x) - log(target_hr_death),
  interval = c(0, 2), tol = 1e-8
)$root

death_trt <- exp(mu_ctrl[2] + delta_d + z_trt[, 2])

recurrence_log_hr_for_shift <- function(delta_r) {
  recurrence_trt <- exp(mu_ctrl[1] + delta_r + z_trt[, 1])
  recurrence_status_trt <- as.integer(recurrence_trt <= death_trt)
  recurrence_followup_trt <- pmin(recurrence_trt, death_trt)
  cox_log_hr(
    recurrence_followup_ctrl, recurrence_status_ctrl,
    recurrence_followup_trt, recurrence_status_trt
  )
}

delta_r <- uniroot(
  function(x) recurrence_log_hr_for_shift(x) - log(target_hr_recurrence),
  interval = c(0, 2.5), tol = 1e-8
)$root

mu_trt <- mu_ctrl + c(delta_r, delta_d)

hr_death <- exp(death_log_hr_for_shift(delta_d))
hr_recurrence <- exp(recurrence_log_hr_for_shift(delta_r))

cat("Standard H1 calibration with Sigma.trt = Sigma.ctrl\n")
cat("Monte Carlo N per arm:", n_per_arm, "\n")
cat(sprintf("delta_R = %.8f; delta_D = %.8f\n", delta_r, delta_d))
cat(sprintf("mu.ctrl = c(%.8f, %.8f)\n", mu_ctrl[1], mu_ctrl[2]))
cat(sprintf("mu.trt  = c(%.8f, %.8f)\n", mu_trt[1], mu_trt[2]))
cat(sprintf(
  "mu.trt  = c(log(%.4f), log(%.4f))\n",
  exp(mu_trt[1]), exp(mu_trt[2])
))
cat(sprintf("recurrence HR = %.6f; death HR = %.6f\n",
            hr_recurrence, hr_death))


standard <- data.frame(arm=c("control","treatment"),mu_log_R=c(mu_ctrl[1],mu_trt[1]),mu_log_D=c(mu_ctrl[2],mu_trt[2]),
  sigma_log_R=sqrt(sigma_common[1,1]),sigma_log_D=sqrt(sigma_common[2,2]),cov_log_R_D=sigma_common[1,2])
write.csv(standard,file.path(output_dir,"standard_h1_parameters.csv"),row.names=FALSE)
