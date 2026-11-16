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

suppressPackageStartupMessages(library(mvtnorm))


days_per_month <- 365.25 / 12
fit_limit <- 84
grid <- seq(0.5, fit_limit, by = 0.5)

baseline <- list(
  Obs = c(mu_r = 3.81071002, sigma_r = 1.820719,
          mu_d = 4.28961452, sigma_d = 1.372863, rho = 0.908254),
  `Lev+5FU` = c(mu_r = 4.74055746, sigma_r = 2.219332,
                mu_d = 4.93960697, sigma_d = 1.864197, rho = 0.945341)
)

step_to_grid <- function(event_times, values, initial, grid) {
  approx(
    x = c(0, event_times),
    y = c(initial, values),
    xout = grid,
    method = "constant",
    f = 0,
    rule = 2,
    ties = "ordered"
  )$y
}

km_curve <- function(time, status, grid) {
  unique_times <- sort(unique(time))
  surv <- 1
  values <- numeric(length(unique_times))
  for (k in seq_along(unique_times)) {
    t <- unique_times[k]
    n_risk <- sum(time >= t)
    d <- sum(time == t & status == 1L)
    surv <- surv * (1 - d / n_risk)
    values[k] <- surv
  }
  step_to_grid(unique_times, values, 1, grid)
}

recurrence_cif_curve <- function(time, state, grid) {
  unique_times <- sort(unique(time))
  event_free <- 1
  cif <- 0
  values <- numeric(length(unique_times))
  for (k in seq_along(unique_times)) {
    t <- unique_times[k]
    n_risk <- sum(time >= t)
    d_r <- sum(time == t & state == "recurrence")
    d_d <- sum(time == t & state == "death")
    cif <- cif + event_free * d_r / n_risk
    event_free <- event_free * (1 - (d_r + d_d) / n_risk)
    values[k] <- cif
  }
  step_to_grid(unique_times, values, 0, grid)
}

empirical_targets <- function(d, grid) {
  r_time <- d$recurrence_time_days / days_per_month
  d_time <- d$death_time_days / days_per_month
  r_status <- as.integer(d$recurrence_status)
  d_status <- as.integer(d$death_status)
  first_time <- pmin(r_time, d_time)
  first_state <- ifelse(
    r_status == 1L & r_time < d_time,
    "recurrence",
    ifelse(d_status == 1L & d_time <= r_time, "death", "censor")
  )
  list(
    death = km_curve(d_time, d_status, grid),
    recurrence = recurrence_cif_curve(first_time, first_state, grid),
    death_weight = vapply(grid, function(t) mean(d_time >= t), numeric(1)),
    recurrence_weight = vapply(
      grid, function(t) mean(first_time >= t), numeric(1)
    )
  )
}

model_death_survival <- function(t, mu_d, sigma_d) {
  plnorm(t, meanlog = mu_d, sdlog = sigma_d, lower.tail = FALSE)
}

model_recurrence_cif <- function(t, mu_r, sigma_r, mu_d, sigma_d, rho) {
  var_v <- sigma_d^2 + sigma_r^2 - 2 * rho * sigma_r * sigma_d
  sd_v <- sqrt(pmax(var_v, 1e-10))
  corr_uv <- (rho * sigma_r * sigma_d - sigma_r^2) /
    (sigma_r * sd_v)
  corr_uv <- max(min(corr_uv, 0.9999), -0.9999)
  z_v <- -(mu_d - mu_r) / sd_v

  vapply(t, function(tt) {
    z_u <- (log(tt) - mu_r) / sigma_r
    joint_lower <- as.numeric(pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(z_u, z_v),
      mean = c(0, 0),
      corr = matrix(c(1, corr_uv, corr_uv, 1), 2, 2),
      algorithm = TVPACK()
    ))
    pmax(pnorm(z_u) - joint_lower, 0)
  }, numeric(1))
}

weighted_mse <- function(delta, weight) {
  sum(weight * delta^2) / sum(weight)
}

calibrate_arm <- function(d, baseline_par, grid) {
  target <- empirical_targets(d, grid)

  death_objective <- function(theta) {
    mu_d <- theta[1]
    sigma_d <- exp(theta[2])
    model <- model_death_survival(grid, mu_d, sigma_d)
    delta <- model - target$death
    weighted_mse(delta, target$death_weight) +
      2 * mean(delta[grid %in% c(12, 24, 36, 60, 84)]^2)
  }

  death_fit <- optim(
    par = c(baseline_par["mu_d"], log(baseline_par["sigma_d"])),
    fn = death_objective,
    method = "L-BFGS-B",
    lower = c(log(5), log(0.25)),
    upper = c(log(500), log(4)),
    control = list(maxit = 2000, factr = 1e5)
  )
  mu_d <- unname(death_fit$par[1])
  sigma_d <- unname(exp(death_fit$par[2]))

  # The paired-data likelihood estimated rho is retained. The recurrence
  # curve alone cannot identify latent recurrence-death dependence.
  rho <- unname(baseline_par["rho"])
  recurrence_objective <- function(theta) {
    mu_r <- theta[1]
    sigma_r <- exp(theta[2])
    model <- model_recurrence_cif(
      grid, mu_r, sigma_r, mu_d, sigma_d, rho
    )
    delta <- model - target$recurrence
    weighted_mse(delta, target$recurrence_weight) +
      2 * mean(delta[grid %in% c(12, 24, 36, 60, 84)]^2) +
      0.75 * max(delta^2) +
      0.0001 * (
        ((mu_r - baseline_par["mu_r"]) / 1.5)^2 +
          (log(sigma_r / baseline_par["sigma_r"]))^2
      )
  }

  starts <- rbind(
    c(baseline_par["mu_r"], log(baseline_par["sigma_r"])),
    c(baseline_par["mu_r"] - 0.25, log(baseline_par["sigma_r"])),
    c(baseline_par["mu_r"] + 0.25, log(baseline_par["sigma_r"])),
    c(baseline_par["mu_r"], log(baseline_par["sigma_r"] * 0.8)),
    c(baseline_par["mu_r"], log(baseline_par["sigma_r"] * 1.2))
  )
  recurrence_fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(
      par = starts[i, ],
      fn = recurrence_objective,
      method = "L-BFGS-B",
      lower = c(log(2), log(0.25)),
      upper = c(log(1000), log(4)),
      control = list(maxit = 3000, factr = 1e5)
    )
  })
  recurrence_fit <- recurrence_fits[[which.min(vapply(
    recurrence_fits, function(x) x$value, numeric(1)
  ))]]
  mu_r <- unname(recurrence_fit$par[1])
  sigma_r <- unname(exp(recurrence_fit$par[2]))

  death_model <- model_death_survival(grid, mu_d, sigma_d)
  recurrence_model <- model_recurrence_cif(
    grid, mu_r, sigma_r, mu_d, sigma_d, rho
  )

  list(
    parameter = c(
      mu_r = mu_r, sigma_r = sigma_r,
      mu_d = mu_d, sigma_d = sigma_d, rho = rho
    ),
    target = target,
    death_model = death_model,
    recurrence_model = recurrence_model,
    diagnostics = data.frame(
      endpoint = c("recurrence", "death"),
      max_absolute_difference = c(
        max(abs(recurrence_model - target$recurrence)),
        max(abs(death_model - target$death))
      ),
      weighted_rmse = c(
        sqrt(weighted_mse(
          recurrence_model - target$recurrence,
          target$recurrence_weight
        )),
        sqrt(weighted_mse(
          death_model - target$death,
          target$death_weight
        ))
      )
    )
  )
}

raw <- colon_patient_pairs(input_file)
raw <- subset(raw, rx %in% c("Obs", "Lev+5FU"))

fits <- list(
  Obs = calibrate_arm(subset(raw, rx == "Obs"), baseline$Obs, grid),
  `Lev+5FU` = calibrate_arm(
    subset(raw, rx == "Lev+5FU"), baseline$`Lev+5FU`, grid
  )
)

parameter_rows <- list()
diagnostic_rows <- list()
for (arm in names(fits)) {
  p <- fits[[arm]]$parameter
  covariance <- p["rho"] * p["sigma_r"] * p["sigma_d"]
  parameter_rows[[arm]] <- data.frame(
    arm = arm,
    mu_R = p["mu_r"],
    sigma_R = p["sigma_r"],
    mu_D = p["mu_d"],
    sigma_D = p["sigma_d"],
    rho = p["rho"],
    var_R = p["sigma_r"]^2,
    covariance_RD = covariance,
    var_D = p["sigma_d"]^2,
    median_R_months = exp(p["mu_r"]),
    median_D_months = exp(p["mu_d"])
  )
  tmp <- fits[[arm]]$diagnostics
  tmp$arm <- arm
  diagnostic_rows[[arm]] <- tmp[, c(
    "arm", "endpoint", "max_absolute_difference", "weighted_rmse"
  )]
}

parameters <- do.call(rbind, parameter_rows)
diagnostics <- do.call(rbind, diagnostic_rows)
write.csv(
  parameters,
  file.path(output_dir, "curve_calibrated_parameters.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(output_dir, "curve_calibrated_theoretical_metrics.csv"),
  row.names = FALSE
)
saveRDS(
  list(parameters = parameters, diagnostics = diagnostics, fits = fits, grid = grid),
  file.path(output_dir, "curve_calibrated_fit.rds")
)

cat("Curve-calibrated bivariate lognormal parameters\n")
print(parameters, row.names = FALSE)
cat("\nTheoretical curve diagnostics through", fit_limit, "months\n")
print(diagnostics, row.names = FALSE)

for (arm in names(fits)) {
  p <- fits[[arm]]$parameter
  cov_rd <- p["rho"] * p["sigma_r"] * p["sigma_d"]
  cat("\n", arm, "\n", sep = "")
  cat(sprintf(
    "mu <- c(log(%.4f), log(%.4f))\n",
    exp(p["mu_r"]), exp(p["mu_d"])
  ))
  cat(sprintf(
    "Sigma <- matrix(c(%.8f, %.8f, %.8f, %.8f), nrow=2, byrow=TRUE)\n",
    p["sigma_r"]^2, cov_rd, cov_rd, p["sigma_d"]^2
  ))
}

# The formal generators use the original source's printed precision:
# medians to 4 decimal places and the covariance matrix to 8 decimal places.
locked <- parameters
locked$median_R_months <- round(locked$median_R_months,4)
locked$median_D_months <- round(locked$median_D_months,4)
locked$mu_R <- log(locked$median_R_months)
locked$mu_D <- log(locked$median_D_months)
locked$var_R <- round(locked$var_R,8)
locked$var_D <- round(locked$var_D,8)
locked$covariance_RD <- round(locked$covariance_RD,8)
locked$sigma_R <- sqrt(locked$var_R)
locked$sigma_D <- sqrt(locked$var_D)
locked$rho <- locked$covariance_RD/(locked$sigma_R*locked$sigma_D)
write.csv(locked,file.path(output_dir,"locked_generator_parameters.csv"),row.names=FALSE)
