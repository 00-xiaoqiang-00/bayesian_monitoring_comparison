#!/usr/bin/env Rscript
# Adopted baseline fit. Numerical functions are retained from the original fit.
# The marginal pseudo-IPD are separate reconstructions; never pair endpoints.
suppressPackageStartupMessages(library(survival))

clean_ipd <- function(ipd) {
  df <- data.frame(time = as.numeric(ipd$time), status = as.integer(ipd$status))
  df <- df[is.finite(df$time) & df$time > 0, , drop = FALSE]
  df
}

km_objects <- function(time, status) {
  fit <- survfit(Surv(time, status) ~ 1, data = data.frame(time = time, status = status))
  med <- as.numeric(summary(fit)$table["median"])
  list(fit = fit, median = med,
       time = c(0, fit$time), surv = c(1, fit$surv),
       n_events = sum(status == 1), n = length(time))
}

km_surv_at <- function(km, t) {
  approx(km$time, km$surv, xout = t, method = "constant", f = 0,
         yleft = 1, yright = min(km$surv, na.rm = TRUE), rule = 2)$y
}

os_milestones <- function(mu, sigma, km, months) {
  data.frame(
    month = months,
    S_model = plnorm(months, meanlog = mu, sdlog = sigma, lower.tail = FALSE),
    S_km = vapply(months, function(t) km_surv_at(km, t), numeric(1)),
    stringsAsFactors = FALSE
  )
}

fit_sigma_fixed_mu <- function(time, status, mu) {
  y <- log(pmax(time, .Machine$double.eps))
  d <- as.integer(status == 1)
  nll <- function(log_s) {
    s <- exp(log_s)
    z <- (y - mu) / s
    ll <- d * (dnorm(z, log = TRUE) - log(s) - y) +
      (1 - d) * pnorm(z, lower.tail = FALSE, log.p = TRUE)
    -sum(ll)
  }
  opt <- optim(0, nll, method = "L-BFGS-B", lower = log(0.05), upper = log(4))
  list(sigma = exp(opt$par), nll = opt$value, conv = opt$convergence)
}

fit_lognormal_mle <- function(time, status) {
  y <- log(pmax(time, .Machine$double.eps))
  d <- as.integer(status == 1)
  nll <- function(par) {
    mu <- par[1]; s <- exp(par[2])
    z <- (y - mu) / s
    ll <- d * (dnorm(z, log = TRUE) - log(s) - y) +
      (1 - d) * pnorm(z, lower.tail = FALSE, log.p = TRUE)
    -sum(ll)
  }
  init <- c(mean(y), log(sd(y) + 0.1))
  opt <- optim(init, nll, method = "L-BFGS-B",
               lower = c(min(y) - 2, log(0.05)),
               upper = c(max(y) + 2, log(4)))
  list(mu = opt$par[1], sigma = exp(opt$par[2]), nll = opt$value)
}

fit_sigma_wls_km <- function(km_time, km_surv, mu, t_max = NULL) {
  ok <- km_time > 0 & km_surv > 0 & km_surv < 1
  if (!is.null(t_max)) ok <- ok & km_time <= t_max
  tt <- km_time[ok]; ss <- km_surv[ok]
  if (length(tt) < 5) return(NA_real_)
  w <- pmax(abs(diff(c(1, ss))), 1e-4); w <- w / sum(w)
  obj <- function(log_s) {
    s <- exp(log_s)
    sum(w * (plnorm(tt, mu, s, lower.tail = FALSE) - ss)^2)
  }
  opt <- optim(0, obj, method = "L-BFGS-B", lower = log(0.05), upper = log(4))
  exp(opt$par)
}

fit_sigma_minimax_km <- function(km_time, km_surv, mu, t_max = NULL) {
  ok <- km_time > 0
  if (!is.null(t_max)) ok <- ok & km_time <= t_max
  tt <- km_time[ok]; ss <- km_surv[ok]
  if (length(tt) < 5) return(NA_real_)
  obj <- function(log_s) {
    max(abs(plnorm(tt, meanlog = mu, sdlog = exp(log_s), lower.tail = FALSE) - ss))
  }
  opt <- optim(0, obj, method = "L-BFGS-B", lower = log(0.05), upper = log(4))
  exp(opt$par)
}

fit_mu_sigma_minimax_km <- function(km_time, km_surv, med_anchor, t_max, med_tol) {
  ok <- km_time > 0
  if (!is.null(t_max)) ok <- ok & km_time <= t_max
  tt <- km_time[ok]; ss <- km_surv[ok]
  if (length(tt) < 5) return(NULL)
  mu_lo <- log(med_anchor * (1 - med_tol))
  mu_hi <- log(med_anchor * (1 + med_tol))
  obj <- function(par) {
    mu <- par[1]; s <- exp(par[2])
    max(abs(plnorm(tt, mu, s, lower.tail = FALSE) - ss))
  }
  opt <- optim(c(log(med_anchor), 0), obj, method = "L-BFGS-B",
               lower = c(mu_lo, log(0.05)), upper = c(mu_hi, log(4)))
  list(mu = opt$par[1], sigma = exp(opt$par[2]), max_abs = opt$value)
}

eval_os_fit <- function(mu, sigma, km, t_main_max, os_max_abs_tol, os_milestone_months) {
  tt <- km$time; ss <- km$surv
  ok <- tt > 0 & tt <= t_main_max
  if (!any(ok)) ok <- tt > 0
  t_use <- tt[ok]; s_km <- ss[ok]
  s_mod <- plnorm(t_use, meanlog = mu, sdlog = sigma, lower.tail = FALSE)
  max_abs <- max(abs(s_mod - s_km))
  integ <- if (length(t_use) >= 2) sum(abs(s_mod - s_km)[-1] * diff(t_use)) else NA_real_
  med_mod <- exp(mu)
  ms <- os_milestones(mu, sigma, km, os_milestone_months)
  list(
    max_abs = max_abs, integ_abs = integ,
    med_model = med_mod, med_km = km$median,
    med_abs_err = abs(med_mod - km$median),
    milestones = ms,
    max_abs_milestone = max(abs(ms$S_model - ms$S_km)),
    pass = is.finite(max_abs) && max_abs <= os_max_abs_tol
  )
}

lock_os_arm <- function(ipd_os, med_anchor, med_lit, label,
                        os_max_abs_tol, os_milestone_months, mu_d_med_tol,
                        t_main_frac = 0.9) {
  df <- clean_ipd(ipd_os)
  km <- km_objects(df$time, df$status)
  t_main <- quantile(df$time, t_main_frac, names = FALSE)
  if (!is.finite(med_anchor) || med_anchor <= 0) med_anchor <- km$median

  candidates <- list()

  # A. μ 钉死锚定中位
  mu0 <- log(med_anchor)
  for (alt in list(
    list(fn = function(...) fit_sigma_fixed_mu(df$time, df$status, mu0)$sigma, name = "fixed_mu_MLE"),
    list(fn = function(...) fit_sigma_wls_km(km$time, km$surv, mu0, t_main), name = "fixed_mu_WLS_KM"),
    list(fn = function(...) fit_sigma_minimax_km(km$time, km$surv, mu0, t_main), name = "fixed_mu_minimax_KM")
  )) {
    sig <- alt$fn()
    if (!is.finite(sig)) next
    candidates[[length(candidates) + 1L]] <- list(
      mu = mu0, sigma = sig, method = alt$name,
      eval = eval_os_fit(mu0, sig, km, t_main, os_max_abs_tol, os_milestone_months)
    )
  }

  # B. μ 在 ±mu_d_med_tol 中位窗内与 σ 联合 minimax
  js <- fit_mu_sigma_minimax_km(km$time, km$surv, med_anchor, t_main, mu_d_med_tol)
  if (!is.null(js)) {
    candidates[[length(candidates) + 1L]] <- list(
      mu = js$mu, sigma = js$sigma,
      method = sprintf("flex_mu_minimax_KM(±%.0f%%)", 100 * mu_d_med_tol),
      eval = eval_os_fit(js$mu, js$sigma, km, t_main, os_max_abs_tol, os_milestone_months)
    )
  }

  # C. 无约束 MLE，再投影到 ±mu_d_med_tol 窗（若超出则夹回边界并重估 σ）
  mle <- fit_lognormal_mle(df$time, df$status)
  med_mle <- exp(mle$mu)
  if (abs(med_mle / med_anchor - 1) <= mu_d_med_tol) {
    candidates[[length(candidates) + 1L]] <- list(
      mu = mle$mu, sigma = mle$sigma, method = "free_MLE_in_tol",
      eval = eval_os_fit(mle$mu, mle$sigma, km, t_main, os_max_abs_tol, os_milestone_months)
    )
  } else {
    mu_clip <- log(med_anchor * (1 + sign(med_mle / med_anchor - 1) * mu_d_med_tol))
    sig_clip <- fit_sigma_fixed_mu(df$time, df$status, mu_clip)$sigma
    candidates[[length(candidates) + 1L]] <- list(
      mu = mu_clip, sigma = sig_clip,
      method = sprintf("MLE_clipped_±%.0f%%", 100 * mu_d_med_tol),
      eval = eval_os_fit(mu_clip, sig_clip, km, t_main, os_max_abs_tol, os_milestone_months)
    )
  }

  passed <- which(vapply(candidates, function(c) isTRUE(c$eval$pass), logical(1)))
  if (length(passed) > 0) {
    i <- passed[which.min(vapply(candidates[passed], function(c) c$eval$max_abs, numeric(1)))]
    best_effort <- FALSE
  } else {
    i <- which.min(vapply(candidates, function(c) c$eval$max_abs, numeric(1)))
    best_effort <- TRUE
  }
  best <- candidates[[i]]

  list(
    label = label, mu = best$mu, sigma = best$sigma, method = best$method,
    med_anchor = med_anchor, med_lit = med_lit, med_km = km$median,
    km = km, t_main = t_main, eval = best$eval,
    best_effort = best_effort,
    n = nrow(df), n_events = sum(df$status == 1)
  )
}
.file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(.file_arg) != 1L) stop("Run this entry with Rscript.")
.script_dir <- dirname(normalizePath(sub("^--file=", "", .file_arg)))
.package_dir <- normalizePath(file.path(.script_dir, "../../.."))
.args <- commandArgs(trailingOnly = TRUE)

# Explicit output destination; resolve existing ancestors and symbolic links
# before checking that the requested directory is outside the submission/code package.
.resolve_output_dir <- function(args, package_dir) {
  if (length(args) == 2L && args[1] == "--output-dir") {
    requested <- args[2]
  } else if (length(args) == 1L && !startsWith(args[1], "--")) {
    requested <- args[1]
  } else {
    stop("Supply --output-dir /absolute/external/directory (or one positional directory).")
  }
  if (!nzchar(requested)) stop("The output directory must not be empty.")
  requested <- path.expand(requested)
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", requested)
  if (!absolute) requested <- file.path(getwd(), requested)
  ancestor <- requested
  suffix <- character()
  while (!file.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) stop("Cannot resolve the output directory.")
    suffix <- c(basename(ancestor), suffix)
    ancestor <- parent
  }
  if (!dir.exists(ancestor)) stop("An output path ancestor is not a directory.")
  resolved <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
  for (part in suffix) {
    if (part == ".") next
    resolved <- if (part == "..") dirname(resolved) else file.path(resolved, part)
    if (file.exists(resolved)) {
      if (!dir.exists(resolved)) stop("An output path ancestor is not a directory.")
      resolved <- normalizePath(resolved, winslash = "/", mustWork = TRUE)
    }
  }
  root <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(dirname(root), "04_Data_and_Results"))) root <- dirname(root)
  output_key <- gsub("\\\\", "/", resolved)
  root_key <- gsub("\\\\", "/", root)
  if (.Platform$OS.type == "windows") {
    output_key <- tolower(output_key)
    root_key <- tolower(root_key)
  }
  output_key <- sub("/+$", "", output_key)
  root_key <- sub("/+$", "", root_key)
  if (identical(output_key, root_key) || startsWith(output_key, paste0(root_key, "/"))) {
    stop("The output directory must be outside the submission/code package.")
  }
  resolved
}
.output_dir <- .resolve_output_dir(.args, .package_dir)

dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(file.path(.package_dir, "inputs/reconstructed_data/fitting_anchors/oak_os_margin/OAK_os_marginal_ipd.rds"))
ipd <- obj$ipd
os <- list()
for (arm in c("trt", "ctrl")) {
  x <- ipd[[paste0("os_", arm)]]$IPD
  d <- clean_ipd(x)
  anchor <- km_objects(d$time, d$status)$median
  med_lit <- if (arm == "trt") 13.8 else 9.6
  os[[arm]] <- lock_os_arm(x, anchor, med_lit, arm,
                          os_max_abs_tol = 0.03,
                          os_milestone_months = c(6, 12, 18),
                          mu_d_med_tol = 0.02)
}
saveRDS(os, file.path(.output_dir, "fitted_os_margin.rds"))
rows <- do.call(rbind, lapply(names(os), function(arm) {
  x <- os[[arm]]
  data.frame(arm = if (arm == "trt") "treatment" else "control",
             mu_D = x$mu, median_D = exp(x$mu), sigma_D = x$sigma,
             adopted_mu_D = log(as.numeric(sprintf("%.4f", exp(x$mu)))),
             adopted_sigma_D = as.numeric(sprintf("%.4f", x$sigma)),
             selected_method = x$method, max_abs_error = x$eval$max_abs,
             passed_original_tolerance = x$eval$pass,
             best_effort = x$best_effort)
}))
write.csv(rows, file.path(.output_dir, "fitted_os_margin.csv"), row.names = FALSE)
print(rows)
