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

suppressPackageStartupMessages({
  library(mvtnorm)
  library(numDeriv)
  library(survival)
})


days_per_month <- 365.25 / 12
analysis_arms <- c("Obs", "Lev+5FU")
landmarks <- c(12, 24, 36, 60, 84)

raw <- colon_patient_pairs(input_file)
dat <- subset(raw, rx %in% analysis_arms)
dat$group <- factor(dat$rx, levels = analysis_arms)
dat$r_time <- dat$recurrence_time_days / days_per_month
dat$d_time <- dat$death_time_days / days_per_month
dat$r_status <- as.integer(dat$recurrence_status)
dat$d_status <- as.integer(dat$death_status)

stopifnot(
  all(dat$r_time > 0),
  all(dat$d_time > 0),
  all(dat$r_status %in% 0:1),
  all(dat$d_status %in% 0:1),
  all(dat$r_time <= dat$d_time)
)

unpack_theta <- function(theta) {
  list(
    mu_r = theta[1],
    sigma_r = exp(theta[2]),
    mu_d = theta[3],
    sigma_d = exp(theta[4]),
    rho = tanh(theta[5])
  )
}

joint_survival <- function(z_r, z_d, rho) {
  vapply(seq_along(z_r), function(i) {
    as.numeric(pmvnorm(
      lower = c(z_r[i], z_d[i]),
      upper = c(Inf, Inf),
      mean = c(0, 0),
      corr = matrix(c(1, rho, rho, 1), 2, 2),
      algorithm = TVPACK()
    ))
  }, numeric(1))
}

neg_loglik <- function(theta, d) {
  p <- unpack_theta(theta)
  if (!is.finite(p$sigma_r) || !is.finite(p$sigma_d) ||
      p$sigma_r < 0.03 || p$sigma_r > 5 ||
      p$sigma_d < 0.03 || p$sigma_d > 5 ||
      abs(p$rho) > 0.995) {
    return(1e12)
  }

  y_r <- log(d$r_time)
  y_d <- log(d$d_time)
  ll <- numeric(nrow(d))

  idx11 <- d$r_status == 1 & d$d_status == 1
  if (any(idx11)) {
    sigma <- matrix(c(
      p$sigma_r^2, p$rho * p$sigma_r * p$sigma_d,
      p$rho * p$sigma_r * p$sigma_d, p$sigma_d^2
    ), 2, 2, byrow = TRUE)
    ll[idx11] <- dmvnorm(
      cbind(y_r[idx11], y_d[idx11]),
      mean = c(p$mu_r, p$mu_d), sigma = sigma, log = TRUE
    ) - y_r[idx11] - y_d[idx11]
  }

  idx10 <- d$r_status == 1 & d$d_status == 0
  if (any(idx10)) {
    cond_mu_d <- p$mu_d + p$rho * p$sigma_d / p$sigma_r *
      (y_r[idx10] - p$mu_r)
    cond_sd_d <- p$sigma_d * sqrt(1 - p$rho^2)
    ll[idx10] <- dnorm(y_r[idx10], p$mu_r, p$sigma_r, log = TRUE) -
      y_r[idx10] +
      pnorm(y_d[idx10], cond_mu_d, cond_sd_d,
            lower.tail = FALSE, log.p = TRUE)
  }

  idx01 <- d$r_status == 0 & d$d_status == 1
  if (any(idx01)) {
    cond_mu_r <- p$mu_r + p$rho * p$sigma_r / p$sigma_d *
      (y_d[idx01] - p$mu_d)
    cond_sd_r <- p$sigma_r * sqrt(1 - p$rho^2)
    ll[idx01] <- dnorm(y_d[idx01], p$mu_d, p$sigma_d, log = TRUE) -
      y_d[idx01] +
      pnorm(y_d[idx01], cond_mu_r, cond_sd_r,
            lower.tail = FALSE, log.p = TRUE)
  }

  idx00 <- d$r_status == 0 & d$d_status == 0
  if (any(idx00)) {
    z_r <- (y_r[idx00] - p$mu_r) / p$sigma_r
    z_d <- (y_d[idx00] - p$mu_d) / p$sigma_d
    probs <- joint_survival(z_r, z_d, p$rho)
    ll[idx00] <- log(pmax(probs, .Machine$double.xmin))
  }

  if (any(!is.finite(ll))) return(1e12)
  -sum(ll)
}

initial_theta <- function(d, rho_start) {
  fit_r <- survreg(Surv(r_time, r_status) ~ 1, data = d,
                   dist = "lognormal")
  fit_d <- survreg(Surv(d_time, d_status) ~ 1, data = d,
                   dist = "lognormal")
  c(
    unname(coef(fit_r)[1]), log(fit_r$scale),
    unname(coef(fit_d)[1]), log(fit_d$scale),
    atanh(rho_start)
  )
}

fit_one_arm <- function(d) {
  rho_starts <- c(-0.25, 0, 0.25, 0.5, 0.75)
  fits <- lapply(rho_starts, function(rho0) {
    start <- initial_theta(d, rho0)
    optim(
      par = start,
      fn = neg_loglik,
      d = d,
      method = "L-BFGS-B",
      lower = c(log(0.05), log(0.05), log(0.05), log(0.05), atanh(-0.98)),
      upper = c(log(500), log(5), log(500), log(5), atanh(0.98)),
      control = list(maxit = 3000, factr = 1e7, pgtol = 1e-8)
    )
  })
  values <- vapply(fits, function(x) x$value, numeric(1))
  best <- fits[[which.min(values)]]

  polished <- optim(
    par = best$par,
    fn = neg_loglik,
    d = d,
    method = "BFGS",
    control = list(maxit = 3000, reltol = 1e-10),
    hessian = TRUE
  )
  if (is.finite(polished$value) && polished$value <= best$value + 1e-5) {
    best <- polished
  }

  h <- if (!is.null(best$hessian)) best$hessian else hessian(neg_loglik, best$par, d = d)
  eig <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  vcov_t <- tryCatch(solve(h), error = function(e) matrix(NA_real_, 5, 5))
  grad_max <- max(abs(grad(neg_loglik, best$par, d = d)))

  list(
    optim = best,
    parameter = unpack_theta(best$par),
    hessian = h,
    hessian_eigen = eig,
    vcov_transformed = vcov_t,
    grad_max = grad_max,
    all_start_values = values
  )
}

delta_se <- function(fit) {
  p <- fit$parameter
  jac <- diag(c(1, p$sigma_r, 1, p$sigma_d, 1 - p$rho^2))
  v <- jac %*% fit$vcov_transformed %*% t(jac)
  sqrt(pmax(diag(v), 0))
}

recurrence_cif_model <- function(t, p) {
  upper <- (log(t) - p$mu_r) / p$sigma_r
  integrand <- function(z_r) {
    cond_mu_d <- p$mu_d + p$rho * p$sigma_d * z_r
    cond_sd_d <- p$sigma_d * sqrt(1 - p$rho^2)
    dnorm(z_r) * pnorm(
      log(exp(p$mu_r + p$sigma_r * z_r)),
      mean = cond_mu_d,
      sd = cond_sd_d,
      lower.tail = FALSE
    )
  }
  integrate(integrand, lower = -8, upper = upper,
            rel.tol = 1e-8, subdivisions = 500L)$value
}

empirical_landmarks <- function(d, times) {
  sf_d <- survfit(Surv(d_time, d_status) ~ 1, data = d)
  os <- summary(sf_d, times = times, extend = TRUE)$surv

  first_time <- ifelse(d$r_status == 1, d$r_time, d$d_time)
  first_event <- ifelse(
    d$r_status == 1, "recurrence",
    ifelse(d$d_status == 1, "death", "censor")
  )
  ms_event <- factor(first_event, levels = c("censor", "recurrence", "death"))
  sf_ms <- survfit(Surv(first_time, ms_event) ~ 1)
  sm <- summary(sf_ms, times = times, extend = TRUE)
  cif_r <- sm$pstate[, "recurrence"]

  data.frame(time_months = times, os_km = os, recurrence_cif_aj = cif_r)
}

groups <- split(dat, dat$group, drop = TRUE)
fits <- lapply(groups, fit_one_arm)

parameter_rows <- list()
diagnostic_rows <- list()
landmark_rows <- list()

for (arm in names(fits)) {
  fit <- fits[[arm]]
  p <- fit$parameter
  se <- delta_se(fit)
  parameter_rows[[arm]] <- data.frame(
    arm = arm,
    parameter = c("mu_R", "sigma_R", "mu_D", "sigma_D", "rho"),
    estimate = c(p$mu_r, p$sigma_r, p$mu_d, p$sigma_d, p$rho),
    se = se,
    lower_95 = c(p$mu_r, p$sigma_r, p$mu_d, p$sigma_d, p$rho) - 1.96 * se,
    upper_95 = c(p$mu_r, p$sigma_r, p$mu_d, p$sigma_d, p$rho) + 1.96 * se
  )

  d <- groups[[arm]]
  diagnostic_rows[[arm]] <- data.frame(
    arm = arm,
    n = nrow(d),
    n_00 = sum(d$r_status == 0 & d$d_status == 0),
    n_01 = sum(d$r_status == 0 & d$d_status == 1),
    n_10 = sum(d$r_status == 1 & d$d_status == 0),
    n_11 = sum(d$r_status == 1 & d$d_status == 1),
    neg_loglik = fit$optim$value,
    aic = 2 * fit$optim$value + 2 * 5,
    convergence = fit$optim$convergence,
    max_abs_gradient = fit$grad_max,
    min_hessian_eigen = min(fit$hessian_eigen),
    median_R_months = exp(p$mu_r),
    median_D_months = exp(p$mu_d),
    prob_R_before_D = recurrence_cif_model(1e6, p)
  )

  lm <- empirical_landmarks(d, landmarks)
  lm$arm <- arm
  lm$os_model <- plnorm(landmarks, p$mu_d, p$sigma_d,
                        lower.tail = FALSE)
  lm$recurrence_cif_model <- vapply(
    landmarks, recurrence_cif_model, numeric(1), p = p
  )
  landmark_rows[[arm]] <- lm[, c(
    "arm", "time_months", "os_km", "os_model",
    "recurrence_cif_aj", "recurrence_cif_model"
  )]

  sigma_matrix <- matrix(c(
    p$sigma_r^2, p$rho * p$sigma_r * p$sigma_d,
    p$rho * p$sigma_r * p$sigma_d, p$sigma_d^2
  ), 2, 2, byrow = TRUE,
  dimnames = list(c("log_R", "log_D"), c("log_R", "log_D")))
  write.csv(sigma_matrix,
            file.path(output_dir, paste0("Sigma_", gsub("\\+", "plus", arm), ".csv")),
            row.names = TRUE)
}

parameters <- do.call(rbind, parameter_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
landmark_check <- do.call(rbind, landmark_rows)

write.csv(parameters, file.path(output_dir, "joint_lognormal_parameters.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "fit_diagnostics.csv"), row.names = FALSE)
write.csv(landmark_check, file.path(output_dir, "landmark_validation.csv"), row.names = FALSE)
saveRDS(fits, file.path(output_dir, "joint_lognormal_fits.rds"))

summary_file <- file.path(output_dir, "fit_summary.txt")
sink(summary_file)
cat("survival_colon bivariate lognormal semi-competing-risk fit\n")
cat("Input:", input_file, "\n")
cat("Time unit: months (365.25/12 days)\n")
cat("Arms: Obs vs Lev+5FU\n\n")
print(diagnostics, row.names = FALSE, digits = 6)
cat("\nParameter estimates\n")
print(parameters, row.names = FALSE, digits = 6)
cat("\nLandmark validation\n")
print(landmark_check, row.names = FALSE, digits = 6)
cat("\nR-ready data-generation parameters\n")
for (arm in names(fits)) {
  p <- fits[[arm]]$parameter
  cov_rd <- p$rho * p$sigma_r * p$sigma_d
  cat("\n", arm, "\n", sep = "")
  cat(sprintf("mu <- c(%.8f, %.8f) # c(log R, log D)\n", p$mu_r, p$mu_d))
  cat(sprintf(
    "Sigma <- matrix(c(%.8f, %.8f, %.8f, %.8f), nrow=2, byrow=TRUE)\n",
    p$sigma_r^2, cov_rd, cov_rd, p$sigma_d^2
  ))
}
sink()

cat(readLines(summary_file), sep = "\n")
