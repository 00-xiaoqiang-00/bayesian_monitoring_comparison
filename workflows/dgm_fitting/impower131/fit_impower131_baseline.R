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

km_weights_at <- function(km, t_pts) {
  fit <- km$fit
  if (is.null(fit) || is.null(fit$std.err) || length(fit$time) < 2) {
    w <- rep(1, length(t_pts))
    return(w / sum(w))
  }
  se <- approx(fit$time, fit$std.err, xout = t_pts, method = "constant",
               f = 0, yleft = fit$std.err[1],
               yright = tail(fit$std.err, 1), rule = 2)$y
  se[!is.finite(se) | se <= 0] <- NA_real_
  w <- 1 / pmax(se, 1e-3)^2
  if (!any(is.finite(w))) w <- rep(1, length(t_pts))
  w[!is.finite(w)] <- min(w[is.finite(w)], na.rm = TRUE)
  w / sum(w)
}

os_milestones <- function(mu, sigma, km, months) {
  data.frame(
    month = months,
    S_model = plnorm(months, meanlog = mu, sdlog = sigma, lower.tail = FALSE),
    S_km = vapply(months, function(t) km_surv_at(km, t), numeric(1)),
    stringsAsFactors = FALSE
  )
}

# 风险集人数：km$fit$n.risk 是右连续阶梯函数（在每个事件/删失时刻之前的人数）
km_nrisk_at <- function(km, t) {
  approx(km$fit$time, km$fit$n.risk, xout = t, method = "constant", f = 0,
         yleft = km$n, yright = 0, rule = 2)$y
}

# ------------------------------------------------------------
# 统一曲线误差评估：目标函数（优化）与最终诊断必须调用同一个函数、
# 同一套 time_grid、同一个 KM 目标，避免"两套不一致的 max|ΔS|"。
# 同时报告：最大误差发生时间、该时刻 KM/DGM 生存率、风险集人数与比例，
# 以及"有效支持区间"（风险集比例 ≥ support_frac_threshold）内的 max|ΔS|。
# ------------------------------------------------------------
evaluate_pfs_curve <- function(mu_R, sig_R, mu_D, sig_D, rho, time_grid, km_target,
                               weights = NULL, support_frac_threshold = 0.10) {
  s_model <- vapply(time_grid, function(t) p_pfs_surv(t, mu_R, sig_R, mu_D, sig_D, rho), numeric(1))
  s_km <- km_surv_at(km_target, time_grid)
  n_risk <- km_nrisk_at(km_target, time_grid)
  frac_risk <- n_risk / km_target$n
  err <- s_model - s_km
  abs_err <- abs(err)

  i_max <- which.max(abs_err)
  ok_support <- frac_risk >= support_frac_threshold
  i_max_supp <- if (any(ok_support)) {
    which(ok_support)[which.max(abs_err[ok_support])]
  } else {
    i_max
  }

  D_km <- if (!is.null(weights)) sum(weights * err^2) else NA_real_

  list(
    time_grid = time_grid, s_model = s_model, s_km = s_km, err = err,
    n_risk = n_risk, frac_risk = frac_risk,
    D_km = D_km,
    max_abs = abs_err[i_max], t_at_max = time_grid[i_max],
    S_km_at_max = s_km[i_max], S_model_at_max = s_model[i_max],
    n_risk_at_max = n_risk[i_max], frac_risk_at_max = frac_risk[i_max],
    max_abs_supported = abs_err[i_max_supp], t_at_max_supported = time_grid[i_max_supp],
    S_km_at_max_supported = s_km[i_max_supp], S_model_at_max_supported = s_model[i_max_supp],
    n_risk_at_max_supported = n_risk[i_max_supp], frac_risk_at_max_supported = frac_risk[i_max_supp],
    support_frac_threshold = support_frac_threshold,
    max_in_thin_tail = frac_risk[i_max] < support_frac_threshold
  )
}

# ------------------------------------------------------------
# 合并（两组共同）删失/失访分布：对重构 OS IPD 做"反向 KM"
# （把删失当事件、死亡当删失），得到潜在删失时间 C_dropout 的生存函数。
# 与我们自己生成的行政删失 C_admin = analysis_time - enroll 独立叠加：
#   C_obs = min(C_admin, C_dropout)
# scale 用于对删失强度做比例校准（S_dropout(t)^scale），scale=1 为原始估计，
# scale>1 删失更快（事件更少），scale<1 删失更慢（事件更多）。
# ------------------------------------------------------------
pooled_dropout_km <- function(ipd_os_trt, ipd_os_ctrl) {
  d1 <- clean_ipd(ipd_os_trt); d2 <- clean_ipd(ipd_os_ctrl)
  df <- rbind(d1, d2)
  fit <- survfit(Surv(time, 1 - status) ~ 1, data = df)
  list(fit = fit, time = c(0, fit$time), surv = c(1, fit$surv),
       median = as.numeric(summary(fit)$table["median"]), n = nrow(df))
}

# u 为预先生成的 Uniform(0,1) 共同随机数（CRN），保证跨 scale/rho 可比
sample_dropout_time <- function(dropout_km, u, scale = 1) {
  ss <- dropout_km$surv^scale
  tt <- dropout_km$time
  ord <- order(tt)
  tt <- tt[ord]; ss <- ss[ord]
  vapply(u, function(ui) {
    idx <- which(ss <= ui)
    if (length(idx) == 0) Inf else tt[idx[1]]
  }, numeric(1))
}

# 用同一套 CRN（z_D、enroll、u_dropout）在给定 dropout scale 下生成 OS 数据，
# 评估：两组事件比例、反向KM中位随访、OS KM曲线拟合、OS HR，联合打分。
calibrate_dropout_scale <- function(cfg, os_trt, os_ctrl, dropout_km, crn,
                                    km_os_t, km_os_c, target_median_fu,
                                    scale_grid = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5),
                                    hr_target = NULL, w_hr = 0.5) {
  if (is.null(hr_target)) hr_target <- cfg$os$HR
  D_t <- exp(os_trt$mu  + os_trt$sigma  * crn$z_D_t)
  D_c <- exp(os_ctrl$mu + os_ctrl$sigma * crn$z_D_c)
  target_evt_t <- cfg$os$events_trt / cfg$n_trt
  target_evt_c <- cfg$os$events_ctrl / cfg$n_ctrl

  eval_scale <- function(scale) {
    C_drop_t <- sample_dropout_time(dropout_km, crn$u_dropout_t, scale)
    C_drop_c <- sample_dropout_time(dropout_km, crn$u_dropout_c, scale)
    C_admin_t <- pmax(cfg$analysis_time - crn$enroll_t, 1e-8)
    C_admin_c <- pmax(cfg$analysis_time - crn$enroll_c, 1e-8)
    C_t <- pmin(C_admin_t, C_drop_t)
    C_c <- pmin(C_admin_c, C_drop_c)
    os_t <- pmin(D_t, C_t); st_t <- as.integer(D_t <= C_t)
    os_c <- pmin(D_c, C_c); st_c <- as.integer(D_c <= C_c)

    evt_t <- mean(st_t); evt_c <- mean(st_c)
    rev_km <- survfit(Surv(c(os_t, os_c), 1 - c(st_t, st_c)) ~ 1)
    med_fu <- as.numeric(summary(rev_km)$table["median"])
    if (!is.finite(med_fu)) med_fu <- max(c(os_t, os_c))

    # 模拟 OS KM vs 真实 IPD OS KM（在真实 IPD 的事件时间上比较）
    km_t <- km_objects(os_t, st_t); km_c <- km_objects(os_c, st_c)
    tt_ipd_t <- km_os_t$time[km_os_t$time > 0]
    tt_ipd_c <- km_os_c$time[km_os_c$time > 0]
    curve_t <- max(abs(km_surv_at(km_t, tt_ipd_t) - km_surv_at(km_os_t, tt_ipd_t)))
    curve_c <- max(abs(km_surv_at(km_c, tt_ipd_c) - km_surv_at(km_os_c, tt_ipd_c)))
    hr <- cox_hr(c(os_t, os_c), c(st_t, st_c), c(rep(1L, length(os_t)), rep(0L, length(os_c))))

    hr_pen <- if (is.finite(hr) && hr > 0) (log(hr) - log(hr_target))^2 else 1
    loss <- (evt_t - target_evt_t)^2 + (evt_c - target_evt_c)^2 +
      ((med_fu - target_median_fu) / target_median_fu)^2 +
      0.5 * (curve_t^2 + curve_c^2) + w_hr * hr_pen
    list(scale = scale, loss = loss, evt_t = evt_t, evt_c = evt_c, med_fu = med_fu,
         curve_t = curve_t, curve_c = curve_c, hr = hr, hr_target = hr_target,
         n_evt_t = sum(st_t), n_evt_c = sum(st_c))
  }

  results <- lapply(scale_grid, eval_scale)
  losses <- vapply(results, `[[`, numeric(1), "loss")
  best <- results[[which.min(losses)]]
  list(best = best, all = results, scale_grid = scale_grid)
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

# 右删失对数正态：联合估 μ,σ（用于 μ 松动搜索的初值）
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

# 固定 μ 附近：对 KM 联合优化 (μ,σ) 最小化 max|ΔS|，约束 |exp(μ)/m_anchor-1|≤med_tol
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

p_pfs_surv <- function(t, mu_R, sig_R, mu_D, sig_D, rho) {
  if (t <= 0) return(1)
  a <- (log(t) - mu_R) / sig_R
  b <- (log(t) - mu_D) / sig_D
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    as.numeric(mvtnorm::pmvnorm(
      lower = c(a, b), upper = c(Inf, Inf),
      mean = c(0, 0),
      corr = matrix(c(1, rho, rho, 1), 2)
    ))
  } else {
    integrand <- function(u) {
      yd <- log(u)
      cond_m <- mu_R + rho * (sig_R / sig_D) * (yd - mu_D)
      cond_s <- sig_R * sqrt(max(1e-12, 1 - rho^2))
      pnorm(log(t), mean = cond_m, sd = cond_s, lower.tail = FALSE) *
        dlnorm(u, meanlog = mu_D, sdlog = sig_D)
    }
    upper <- exp(mu_D + 8 * sig_D)
    tryCatch(
      integrate(integrand, lower = t, upper = upper, rel.tol = 1e-4)$value,
      error = function(e) plnorm(t, mu_R, sig_R, lower.tail = FALSE) *
        plnorm(t, mu_D, sig_D, lower.tail = FALSE)
    )
  }
}

# 硬约束：给定 σ_R，求 μ_R 使 S_PFS(m_target)=0.5
solve_mu_R_for_median <- function(sig_R, mu_D, sig_D, rho, m_target) {
  s_d <- plnorm(m_target, mu_D, sig_D, lower.tail = FALSE)
  if (s_d < 0.5 - 1e-8) {
    return(list(mu_R = Inf, ok = FALSE, S_at_m = s_d))
  }
  f <- function(mu_R) {
    p_pfs_surv(m_target, mu_R, sig_R, mu_D, sig_D, rho) - 0.5
  }
  lo <- log(m_target) - 6
  hi <- log(m_target) + 8
  flo <- f(lo); fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi)) {
    return(list(mu_R = log(m_target), ok = FALSE, S_at_m = NA_real_))
  }
  if (flo > 0) {
    for (k in 1:10) {
      lo <- lo - 2
      flo <- f(lo)
      if (flo <= 0) break
    }
  }
  if (fhi < 0) {
    for (k in 1:10) {
      hi <- hi + 2
      fhi <- f(hi)
      if (fhi >= 0) break
    }
  }
  if (flo > 0 || fhi < 0) {
    return(list(mu_R = log(m_target), ok = FALSE, S_at_m = p_pfs_surv(m_target, log(m_target), sig_R, mu_D, sig_D, rho)))
  }
  rt <- tryCatch(uniroot(f, c(lo, hi), tol = 1e-8), error = function(e) NULL)
  if (is.null(rt)) {
    return(list(mu_R = log(m_target), ok = FALSE, S_at_m = NA_real_))
  }
  list(mu_R = rt$root, ok = TRUE,
       S_at_m = p_pfs_surv(m_target, rt$root, sig_R, mu_D, sig_D, rho))
}

gen_RD <- function(mu_R, sig_R, mu_D, sig_D, rho, z_D, z_E) {
  D <- exp(mu_D + sig_D * z_D)
  R <- exp(mu_R + sig_R * (rho * z_D + sqrt(pmax(0, 1 - rho^2)) * z_E))
  list(R = R, D = D)
}

# 统一删失：C_eff = min(analysis_time - enroll, C_rand)；C_rand = Inf → 仅行政删失
apply_censor <- function(R, D, enroll, C_rand, analysis_time) {
  C_admin <- pmax(analysis_time - enroll, 1e-8)
  C <- pmin(C_admin, C_rand)
  pfs_latent <- pmin(R, D)
  data.frame(
    pfs_time = pmin(pfs_latent, C),
    pfs_status = as.integer(pfs_latent <= C),
    os_time = pmin(D, C),
    os_status = as.integer(D <= C),
    R = R, D = D, C = C, C_admin = C_admin, C_rand = C_rand,
    enroll = enroll
  )
}

cox_hr <- function(time, status, arm) {
  fit <- tryCatch(coxph(Surv(time, status) ~ arm), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(exp(coef(fit)[1]))
}

# 统一试验规模验证：同一套 D / enroll / C_rand → OS（及可选 PFS）指标
# 保证第1–3步与第5步 OS HR / 事件数完全一致
trial_scale_metrics <- function(D_t, D_c, enroll_t, enroll_c, C_rand_t, C_rand_c,
                                analysis_time, R_t = NULL, R_c = NULL) {
  if (is.null(R_t)) R_t <- D_t
  if (is.null(R_c)) R_c <- D_c
  obs_t <- apply_censor(R_t, D_t, enroll_t, C_rand_t, analysis_time)
  obs_c <- apply_censor(R_c, D_c, enroll_c, C_rand_c, analysis_time)
  arm <- c(rep(1L, length(D_t)), rep(0L, length(D_c)))
  list(
    obs_t = obs_t, obs_c = obs_c,
    n_os_evt_t = sum(obs_t$os_status), n_os_evt_c = sum(obs_c$os_status),
    n_pfs_evt_t = sum(obs_t$pfs_status), n_pfs_evt_c = sum(obs_c$pfs_status),
    med_os_t = km_objects(obs_t$os_time, obs_t$os_status)$median,
    med_os_c = km_objects(obs_c$os_time, obs_c$os_status)$median,
    med_pfs_t = km_objects(obs_t$pfs_time, obs_t$pfs_status)$median,
    med_pfs_c = km_objects(obs_c$pfs_time, obs_c$pfs_status)$median,
    HR_os = cox_hr(c(obs_t$os_time, obs_c$os_time),
                   c(obs_t$os_status, obs_c$os_status), arm),
    HR_pfs = cox_hr(c(obs_t$pfs_time, obs_c$pfs_time),
                    c(obs_t$pfs_status, obs_c$pfs_status), arm)
  )
}

# censor_mode: "admin_only" (C_rand=Inf) | "uniform" (Unif(0,cens_upper)) |
#              "km_dropout" (合并两组重构IPD反向KM删失分布，按 dropout_scale 校准强度)
#
# 重要：为保证与旧版逐字节一致的 CRN（不因本次重构而意外改变 admin_only /
# uniform 试验如 IMpower131 的随机数流），仅在 censor_mode=="km_dropout" 时才
# 额外抽取 u_dropout；其余模式下的 RNG 调用序列与重构前完全相同。
`%||%` <- function(a, b) if (is.null(a)) b else a

draw_one_arm_crn <- function(n, enroll_max, mode, cens_upper, dropout_km, dropout_scale) {
  z_D <- rnorm(n); z_E <- rnorm(n)
  enroll <- runif(n, 0, enroll_max)
  if (mode == "km_dropout") {
    u_dropout <- runif(n)
    C_rand <- sample_dropout_time(dropout_km, u_dropout, dropout_scale)
  } else if (mode == "uniform") {
    u_dropout <- NULL
    C_rand <- runif(n, 0, cens_upper)
  } else {
    u_dropout <- NULL
    C_rand <- rep(Inf, n)
  }
  list(z_D = z_D, z_E = z_E, enroll = enroll, u_dropout = u_dropout, C_rand = C_rand)
}

make_crn <- function(cfg, seed, dropout_km = NULL, dropout_scale = 1) {
  mode <- cfg$censor_mode %||% (if (isTRUE(cfg$use_random_censor)) "uniform" else "admin_only")
  set.seed(seed)
  at <- draw_one_arm_crn(cfg$n_trt, cfg$enroll_max, mode, cfg$cens_upper, dropout_km, dropout_scale)
  ac <- draw_one_arm_crn(cfg$n_ctrl, cfg$enroll_max, mode, cfg$cens_upper, dropout_km, dropout_scale)
  list(
    z_D_t = at$z_D, z_E_t = at$z_E, enroll_t = at$enroll,
    u_dropout_t = at$u_dropout, C_rand_t = at$C_rand,
    z_D_c = ac$z_D, z_E_c = ac$z_E, enroll_c = ac$enroll,
    u_dropout_c = ac$u_dropout, C_rand_c = ac$C_rand
  )
}

make_hr_crn <- function(cfg, seed, n = NULL, dropout_km = NULL, dropout_scale = 1) {
  if (is.null(n)) n <- cfg$calibrate$n_hr
  mode <- cfg$censor_mode %||% (if (isTRUE(cfg$use_random_censor)) "uniform" else "admin_only")
  set.seed(seed + 99L)
  at <- draw_one_arm_crn(n, cfg$enroll_max, mode, cfg$cens_upper, dropout_km, dropout_scale)
  ac <- draw_one_arm_crn(n, cfg$enroll_max, mode, cfg$cens_upper, dropout_km, dropout_scale)
  list(
    z_D_t = at$z_D, z_E_t = at$z_E, enroll_t = at$enroll,
    u_dropout_t = at$u_dropout, C_rand_t = at$C_rand,
    z_D_c = ac$z_D, z_E_c = ac$z_E, enroll_c = ac$enroll,
    u_dropout_c = ac$u_dropout, C_rand_c = ac$C_rand
  )
}

fmt_mu <- function(mu) sprintf("log(%.4f)", exp(mu))
fmt_pair <- function(mu) sprintf("c(%s, %s)", fmt_mu(mu[1]), fmt_mu(mu[2]))
make_Sigma <- function(sigR, sigD, rho) {
  matrix(c(sigR^2, rho * sigR * sigD, rho * sigR * sigD, sigD^2), 2, byrow = TRUE)
}

ipd_unadj_hr <- function(ipd_trt, ipd_ctrl) {
  d1 <- clean_ipd(ipd_trt); d0 <- clean_ipd(ipd_ctrl)
  d1$arm <- 1L; d0$arm <- 0L
  df <- rbind(d1, d0)
  cox_hr(df$time, df$status, df$arm)
}

# ------------------------------------------------------------
# 单臂：锁定 OS（允许 μ_D 在锚定中位 ±mu_d_med_tol 内松动）
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 给定锁定 OS + ρ：硬约束中位 + 联合优化两臂 σ_R（含 PFS HR 惩罚）
# 两阶段：① 仅 KM（快）→ ② 加入大样本 HR 惩罚微调
# L_R = w_S*(D_KM,trt + D_KM,ctrl) + w_H*(log HR_DGM - log HR_target)^2
# ------------------------------------------------------------
# 由 KM 事件时间构造"统一时间网格"：不做 t_cap 截断（诊断需要看到全程尾部），
# 仅在点数过多时按密度保序抽稀到 max_n，且这套网格会被目标函数与最终诊断
# 同时复用，从根本上避免"优化用一套网格、诊断用另一套网格"的不一致。
build_common_grid <- function(km, max_n = 80) {
  tt <- sort(unique(km$time))
  tt <- tt[tt > 0]
  if (length(tt) > max_n) {
    tt <- tt[round(seq(1, length(tt), length.out = max_n))]
  }
  tt
}

implied_pfs_median <- function(mu_R, sig_R, mu_D, sig_D, rho, guess) {
  f <- function(logt) p_pfs_surv(exp(logt), mu_R, sig_R, mu_D, sig_D, rho) - 0.5
  lo <- log(guess) - 6; hi <- log(guess) + 8
  flo <- tryCatch(f(lo), error = function(e) NA_real_)
  fhi <- tryCatch(f(hi), error = function(e) NA_real_)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  rt <- tryCatch(uniroot(f, c(lo, hi), tol = 1e-6), error = function(e) NULL)
  if (is.null(rt)) return(NA_real_)
  exp(rt$root)
}

# median_mode = "hard"：给定 σ_R 数值求根使 S_PFS(m_anchor)=0.5（原方法，向后兼容）
# median_mode = "soft"：μ_R,σ_R 均为自由参数，中位数偏差通过软惩罚
#   w_m*(max(0, |rel误差|-tol))^2 约束在 tol（如5%）以内，不强制相等。
# 目标函数与最终诊断统一调用 evaluate_pfs_curve()，使用同一套 build_common_grid()
# 网格，从根本上保证 max|ΔS|_opt 与 max|ΔS|(诊断) 一致。
calibrate_R_given_rho <- function(rho, os_trt, os_ctrl,
                                  ipd_pfs_trt, ipd_pfs_ctrl,
                                  cfg, crn, hr_crn,
                                  med_pfs_anchor_trt, med_pfs_anchor_ctrl,
                                  hr_pfs_target, pfs_milestone_months,
                                  median_mode = "hard",
                                  support_frac_threshold = 0.10) {
  cal <- cfg$calibrate
  w_s <- cal$w_s; w_max <- cal$w_max; w_h <- cal$w_h; w_m <- cal$w_m %||% 0
  med_tol <- cal$pfs_med_rel_tol
  hr_pref_lo <- cal$hr_pref_lo; hr_pref_hi <- cal$hr_pref_hi

  km_pfs_trt  <- km_objects(clean_ipd(ipd_pfs_trt)$time,  clean_ipd(ipd_pfs_trt)$status)
  km_pfs_ctrl <- km_objects(clean_ipd(ipd_pfs_ctrl)$time, clean_ipd(ipd_pfs_ctrl)$status)

  # 统一网格：目标函数（优化）与最终诊断都用这同一套 grid_t / grid_c
  grid_t <- build_common_grid(km_pfs_trt)
  grid_c <- build_common_grid(km_pfs_ctrl)
  w_trt  <- km_weights_at(km_pfs_trt, grid_t)
  w_ctrl <- km_weights_at(km_pfs_ctrl, grid_c)

  curve_eval_t <- function(mu_R, sig_R) {
    evaluate_pfs_curve(mu_R, sig_R, os_trt$mu, os_trt$sigma, rho, grid_t, km_pfs_trt,
                       weights = w_trt, support_frac_threshold = support_frac_threshold)
  }
  curve_eval_c <- function(mu_R, sig_R) {
    evaluate_pfs_curve(mu_R, sig_R, os_ctrl$mu, os_ctrl$sigma, rho, grid_c, km_pfs_ctrl,
                       weights = w_ctrl, support_frac_threshold = support_frac_threshold)
  }

  hr_from_params <- function(mu_Rt, sig_t, mu_Rc, sig_c) {
    rd_t <- gen_RD(mu_Rt, sig_t, os_trt$mu, os_trt$sigma, rho,
                   hr_crn$z_D_t, hr_crn$z_E_t)
    rd_c <- gen_RD(mu_Rc, sig_c, os_ctrl$mu, os_ctrl$sigma, rho,
                   hr_crn$z_D_c, hr_crn$z_E_c)
    trial_scale_metrics(
      rd_t$D, rd_c$D, hr_crn$enroll_t, hr_crn$enroll_c,
      hr_crn$C_rand_t, hr_crn$C_rand_c, cfg$analysis_time,
      rd_t$R, rd_c$R
    )$HR_pfs
  }

  if (median_mode == "hard") {
    curve_parts <- function(sig_t, sig_c) {
      sol_t <- solve_mu_R_for_median(sig_t, os_trt$mu, os_trt$sigma, rho, med_pfs_anchor_trt)
      sol_c <- solve_mu_R_for_median(sig_c, os_ctrl$mu, os_ctrl$sigma, rho, med_pfs_anchor_ctrl)
      if (!isTRUE(sol_t$ok) || !isTRUE(sol_c$ok)) return(NULL)
      ev_t <- curve_eval_t(sol_t$mu_R, sig_t)
      ev_c <- curve_eval_c(sol_c$mu_R, sig_c)
      list(mu_Rt = sol_t$mu_R, sig_t = sig_t, mu_Rc = sol_c$mu_R, sig_c = sig_c,
           ev_t = ev_t, ev_c = ev_c,
           D_km_t = ev_t$D_km, D_km_c = ev_c$D_km,
           max_abs_t = ev_t$max_abs, max_abs_c = ev_c$max_abs,
           med_rel_t = 0, med_rel_c = 0, S_at_m_t = sol_t$S_at_m, S_at_m_c = sol_c$S_at_m)
    }
    med_pen <- function(cp) 0
    npar <- 2L
    par0 <- c(log(os_trt$sigma), log(os_ctrl$sigma))
    lower0 <- c(log(0.15), log(0.15)); upper0 <- c(log(3.5), log(3.5))
    make_par <- function(par) list(sig_t = exp(par[1]), sig_c = exp(par[2]))
    call_curve <- function(par) curve_parts(exp(par[1]), exp(par[2]))
  } else {
    curve_parts <- function(mu_Rt, sig_t, mu_Rc, sig_c) {
      ev_t <- curve_eval_t(mu_Rt, sig_t)
      ev_c <- curve_eval_c(mu_Rc, sig_c)
      med_model_t <- implied_pfs_median(mu_Rt, sig_t, os_trt$mu, os_trt$sigma, rho, med_pfs_anchor_trt)
      med_model_c <- implied_pfs_median(mu_Rc, sig_c, os_ctrl$mu, os_ctrl$sigma, rho, med_pfs_anchor_ctrl)
      med_rel_t <- if (is.finite(med_model_t)) abs(med_model_t / med_pfs_anchor_trt - 1) else 1
      med_rel_c <- if (is.finite(med_model_c)) abs(med_model_c / med_pfs_anchor_ctrl - 1) else 1
      list(mu_Rt = mu_Rt, sig_t = sig_t, mu_Rc = mu_Rc, sig_c = sig_c,
           ev_t = ev_t, ev_c = ev_c,
           D_km_t = ev_t$D_km, D_km_c = ev_c$D_km,
           max_abs_t = ev_t$max_abs, max_abs_c = ev_c$max_abs,
           med_rel_t = med_rel_t, med_rel_c = med_rel_c,
           S_at_m_t = ev_t$s_model[which.min(abs(grid_t - med_pfs_anchor_trt))],
           S_at_m_c = ev_c$s_model[which.min(abs(grid_c - med_pfs_anchor_ctrl))])
    }
    med_pen <- function(cp) {
      w_m * (max(0, cp$med_rel_t - med_tol)^2 + max(0, cp$med_rel_c - med_tol)^2)
    }
    npar <- 4L
    par0 <- c(log(med_pfs_anchor_trt), log(os_trt$sigma),
             log(med_pfs_anchor_ctrl), log(os_ctrl$sigma))
    lower0 <- c(log(med_pfs_anchor_trt) - 1, log(0.15),
               log(med_pfs_anchor_ctrl) - 1, log(0.15))
    upper0 <- c(log(med_pfs_anchor_trt) + 1, log(3.5),
               log(med_pfs_anchor_ctrl) + 1, log(3.5))
    make_par <- function(par) list(mu_Rt = par[1], sig_t = exp(par[2]),
                                   mu_Rc = par[3], sig_c = exp(par[4]))
    call_curve <- function(par) curve_parts(par[1], exp(par[2]), par[3], exp(par[4]))
  }

  # 多目标损失；HR 落在本试验偏好带内轻奖励；软模式下加中位数超容差惩罚
  joint_loss <- function(cp, hr) {
    L <- w_s * (cp$D_km_t + cp$D_km_c) +
      w_max * (cp$max_abs_t + cp$max_abs_c) +
      w_h * (log(hr) - log(hr_pfs_target))^2 +
      med_pen(cp)
    if (is.finite(hr) && hr >= hr_pref_lo && hr <= hr_pref_hi) {
      L <- L - 0.02
    }
    L
  }

  # 阶段1：KM + max|ΔS|（+软模式中位惩罚），不含 HR（更快）
  obj_km <- function(par) {
    cp <- call_curve(par)
    if (is.null(cp)) return(1e6)
    w_s * (cp$D_km_t + cp$D_km_c) + w_max * (cp$max_abs_t + cp$max_abs_c) + med_pen(cp)
  }
  opt1 <- optim(par0, obj_km, method = "L-BFGS-B",
                lower = lower0, upper = upper0,
                control = list(maxit = 60))

  # 阶段2：局部网格 + 完整多目标（含 HR）
  base <- opt1$par
  n_lv <- if (npar <= 2) 7L else 3L
  grid_list <- lapply(seq_len(npar), function(i) {
    unique(pmin(pmax(base[i] + seq(-0.7, 0.7, length.out = n_lv), lower0[i]), upper0[i]))
  })
  grid_df <- expand.grid(grid_list)
  best_par <- base
  best_val <- Inf
  for (i in seq_len(nrow(grid_df))) {
    par_i <- as.numeric(grid_df[i, ])
    cp <- call_curve(par_i)
    if (is.null(cp)) next
    pp <- make_par(par_i)
    hr <- hr_from_params(if (!is.null(pp$mu_Rt)) pp$mu_Rt else cp$mu_Rt, pp$sig_t,
                         if (!is.null(pp$mu_Rc)) pp$mu_Rc else cp$mu_Rc, pp$sig_c)
    if (!is.finite(hr) || hr <= 0) next
    val <- joint_loss(cp, hr)
    if (val < best_val) {
      best_val <- val
      best_par <- par_i
    }
  }
  obj_full <- function(par) {
    cp <- call_curve(par)
    if (is.null(cp)) return(1e6)
    hr <- hr_from_params(cp$mu_Rt, cp$sig_t, cp$mu_Rc, cp$sig_c)
    if (!is.finite(hr) || hr <= 0) return(1e6)
    joint_loss(cp, hr)
  }
  opt2 <- optim(best_par, obj_full, method = "L-BFGS-B",
                lower = pmax(best_par - 0.35, lower0),
                upper = pmin(best_par + 0.35, upper0),
                control = list(maxit = 25))
  cp <- call_curve(opt2$par)
  hr_dgm <- hr_from_params(cp$mu_Rt, cp$sig_t, cp$mu_Rc, cp$sig_c)
  hr_pen <- (log(hr_dgm) - log(hr_pfs_target))^2
  loss <- joint_loss(cp, hr_dgm)

  # 最终诊断：与目标函数使用同一 evaluate_pfs_curve() 结果（ev_t/ev_c），
  # 不再另外重新计算一遍——从根本上消除"两套不一致 max|ΔS|"的可能。
  pfs_diag_from_eval <- function(ev, mu_R, sig_R, os, med_target) {
    ms <- data.frame(
      month = pfs_milestone_months,
      S_model = vapply(pfs_milestone_months, function(t)
        p_pfs_surv(t, mu_R, sig_R, os$mu, os$sigma, rho), numeric(1)),
      S_km = vapply(pfs_milestone_months, function(t) km_surv_at(
        list(time = ev$time_grid, surv = ev$s_km), t), numeric(1))
    )
    ms$abs_err <- abs(ms$S_model - ms$S_km)
    list(
      med_target = med_target,
      S_at_m = p_pfs_surv(med_target, mu_R, sig_R, os$mu, os$sigma, rho),
      max_abs = ev$max_abs, t_at_max = ev$t_at_max,
      S_km_at_max = ev$S_km_at_max, S_model_at_max = ev$S_model_at_max,
      n_risk_at_max = ev$n_risk_at_max, frac_risk_at_max = ev$frac_risk_at_max,
      max_abs_supported = ev$max_abs_supported, t_at_max_supported = ev$t_at_max_supported,
      S_km_at_max_supported = ev$S_km_at_max_supported,
      S_model_at_max_supported = ev$S_model_at_max_supported,
      n_risk_at_max_supported = ev$n_risk_at_max_supported,
      frac_risk_at_max_supported = ev$frac_risk_at_max_supported,
      max_in_thin_tail = ev$max_in_thin_tail,
      max_milestone = max(ms$abs_err),
      milestones = ms
    )
  }
  diag_t <- pfs_diag_from_eval(cp$ev_t, cp$mu_Rt, cp$sig_t, os_trt, med_pfs_anchor_trt)
  diag_c <- pfs_diag_from_eval(cp$ev_c, cp$mu_Rc, cp$sig_c, os_ctrl, med_pfs_anchor_ctrl)

  D_t <- exp(os_trt$mu  + os_trt$sigma  * crn$z_D_t)
  D_c <- exp(os_ctrl$mu + os_ctrl$sigma * crn$z_D_c)
  rd_t <- gen_RD(cp$mu_Rt, cp$sig_t, os_trt$mu, os_trt$sigma, rho,
                 crn$z_D_t, crn$z_E_t)
  rd_c <- gen_RD(cp$mu_Rc, cp$sig_c, os_ctrl$mu, os_ctrl$sigma, rho,
                 crn$z_D_c, crn$z_E_c)
  met <- trial_scale_metrics(
    D_t, D_c, crn$enroll_t, crn$enroll_c, crn$C_rand_t, crn$C_rand_c,
    cfg$analysis_time, rd_t$R, rd_c$R
  )

  med_rel_t <- abs(met$med_pfs_t / med_pfs_anchor_trt - 1)
  med_rel_c <- abs(met$med_pfs_c / med_pfs_anchor_ctrl - 1)

  list(
    rho = rho,
    L_R = loss, D_km_t = cp$D_km_t, D_km_c = cp$D_km_c,
    max_abs_opt_t = cp$max_abs_t, max_abs_opt_c = cp$max_abs_c,
    hr_pen = hr_pen,
    mu_R_trt = cp$mu_Rt, sig_R_trt = cp$sig_t,
    mu_R_ctrl = cp$mu_Rc, sig_R_ctrl = cp$sig_c,
    S_pfs_at_m_trt = cp$S_at_m_t, S_pfs_at_m_ctrl = cp$S_at_m_c,
    median_mode = median_mode,
    med_model_rel_err_t = cp$med_rel_t, med_model_rel_err_c = cp$med_rel_c,
    HR_pfs_dgm = hr_dgm,
    hr_pfs_target = hr_pfs_target,
    hr_pfs_rel_err = abs(hr_dgm / hr_pfs_target - 1),
    hr_in_pref_band = is.finite(hr_dgm) && hr_dgm >= hr_pref_lo && hr_dgm <= hr_pref_hi,
    med_pfs_trt = met$med_pfs_t, med_pfs_ctrl = met$med_pfs_c,
    med_rel_err_t = med_rel_t, med_rel_err_c = med_rel_c,
    med_os_trt = met$med_os_t, med_os_ctrl = met$med_os_c,
    n_os_evt_trt = met$n_os_evt_t, n_os_evt_ctrl = met$n_os_evt_c,
    n_pfs_evt_trt = met$n_pfs_evt_t, n_pfs_evt_ctrl = met$n_pfs_evt_c,
    HR_pfs_trial = met$HR_pfs, HR_os = met$HR_os,
    D_trt_mean = mean(D_t), D_ctrl_mean = mean(D_c),
    pfs_diag_trt = diag_t, pfs_diag_ctrl = diag_c,
    os_fingerprint = c(n_t = met$n_os_evt_t, n_c = met$n_os_evt_c, HR = met$HR_os)
  )
}

run_os_first_calibration <- function(trial_name, cfg,
                                     rho_grid,
                                     anchor_source = "ipd",
                                     out_subdir = NULL,
                                     seed = 2026L,
                                     os_max_abs_tol = 0.03,
                                     os_milestone_months = c(6, 12, 18),
                                     pfs_milestone_months = c(6, 12, 18),
                                     mu_d_med_tol = 0.02,
                                     median_mode = "hard",
                                     support_frac_threshold = 0.10,
                                     rho_main = NULL,
                                     target_median_fu = NULL,
                                     dropout_scale_grid = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5)) {
  if (is.null(out_subdir)) {
    out_subdir <- if (anchor_source == "ipd") "os_first_method_ipd_anchor" else "os_first_method"
  }
  out_dir <- cfg$output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n========== %s (anchor=%s) ==========\n", trial_name, anchor_source))
  obj <- readRDS(cfg$rds)
  ipd <- obj$ipd

  km_os_t  <- km_objects(clean_ipd(ipd$os_trt$IPD)$time,  clean_ipd(ipd$os_trt$IPD)$status)
  km_os_c  <- km_objects(clean_ipd(ipd$os_ctrl$IPD)$time, clean_ipd(ipd$os_ctrl$IPD)$status)
  km_pfs_t <- km_objects(clean_ipd(ipd$pfs_trt$IPD)$time,  clean_ipd(ipd$pfs_trt$IPD)$status)
  km_pfs_c <- km_objects(clean_ipd(ipd$pfs_ctrl$IPD)$time, clean_ipd(ipd$pfs_ctrl$IPD)$status)

  if (anchor_source == "ipd") {
    med_os_anchor_t  <- km_os_t$median
    med_os_anchor_c  <- km_os_c$median
    med_pfs_anchor_t <- km_pfs_t$median
    med_pfs_anchor_c <- km_pfs_c$median
  } else {
    med_os_anchor_t  <- cfg$os$median_trt
    med_os_anchor_c  <- cfg$os$median_ctrl
    med_pfs_anchor_t <- cfg$pfs$median_trt
    med_pfs_anchor_c <- cfg$pfs$median_ctrl
  }

  report <- character()
  add <- function(...) report <<- c(report, paste0(...))

  add("======================================================================")
  add(sprintf("试验：%s | OS先锁 → 硬约束中位校准 R（PFS=min(R,D)）", trial_name))
  add(sprintf("锚定来源：%s", ifelse(anchor_source == "ipd", "IPD重建KM中位", "文献报告中位")))
  add(sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  add(sprintf(
    "删失：analysis_time=%g | enroll~U(0,%g) | C_admin=analysis_time-enroll | 随机失访=%s",
    cfg$analysis_time, cfg$enroll_max,
    ifelse(isTRUE(cfg$use_random_censor),
           sprintf("ON cens_upper=%g", cfg$cens_upper), "OFF (C_rand=Inf)")
  ))
  add("说明：enroll_max 是入组窗口，不是 analysis_time；勿用 runif(0, analysis_time) 当行政删失。")
  add(sprintf(
    "锚定中位 OS  trt/ctrl = %.4f / %.4f（文献 %.1f / %.1f）",
    med_os_anchor_t, med_os_anchor_c, cfg$os$median_trt, cfg$os$median_ctrl
  ))
  add(sprintf(
    "锚定中位 PFS trt/ctrl = %.4f / %.4f（文献 %.1f / %.1f）",
    med_pfs_anchor_t, med_pfs_anchor_c, cfg$pfs$median_trt, cfg$pfs$median_ctrl
  ))
  add("共同随机数：固定 z_D / z_E / enroll / C_rand")
  add("======================================================================\n")

  # ---- 第1–3步：锁 OS ----
  add("--- 第1–3步：校准死亡分布并验证 ---")
  add(sprintf("μ_D 允许相对锚定中位 ±%.0f%% 松动以改善曲线；ρ 网格：%s",
              100 * mu_d_med_tol, paste(sprintf("%.2f", rho_grid), collapse = ", ")))
  os_trt  <- lock_os_arm(ipd$os_trt$IPD,  med_os_anchor_t, cfg$os$median_trt,  "trt",
                        os_max_abs_tol, os_milestone_months, mu_d_med_tol)
  os_ctrl <- lock_os_arm(ipd$os_ctrl$IPD, med_os_anchor_c, cfg$os$median_ctrl, "ctrl",
                        os_max_abs_tol, os_milestone_months, mu_d_med_tol)

  for (os in list(os_trt, os_ctrl)) {
    add(sprintf(
      "[%s] μ_D=%s  σ_D=%.4f  method=%s%s",
      os$label, fmt_mu(os$mu), os$sigma, os$method,
      ifelse(os$best_effort, "  [无合格→最优锁定]", "")
    ))
    add(sprintf(
      "      中位: 模型=%.3f  IPD_KM=%.3f  锚定=%.3f  文献=%.1f  |模型/锚定-1|=%.2f%%",
      os$eval$med_model, os$med_km, os$med_anchor, os$med_lit,
      100 * abs(os$eval$med_model / os$med_anchor - 1)
    ))
    add(sprintf(
      "      主窗 t≤%.1f: max|S_model-S_KM|=%.4f  积分|差|=%.3f  %s (容差≤%.2f)",
      os$t_main, os$eval$max_abs, os$eval$integ_abs,
      ifelse(os$eval$pass, "PASS", "FAIL"), os_max_abs_tol
    ))
    ms <- os$eval$milestones
    for (m in os_milestone_months) {
      add(sprintf(
        "      %2dm OS: 模型=%.3f  KM=%.3f  |Δ|=%.3f",
        m, ms$S_model[ms$month == m], ms$S_km[ms$month == m],
        abs(ms$S_model[ms$month == m] - ms$S_km[ms$month == m])
      ))
    }
    add(sprintf("      IPD OS 事件数=%d / n=%d (%.1f%%)",
                os$n_events, os$n, 100 * os$n_events / os$n))
  }

  # ---- 可选：从重构 OS IPD 估计合并（两组共同）反向KM删失分布，
  # 再按 event%/中位随访/OS KM/OS HR 对删失强度做网格校准 ----
  censor_mode <- cfg$censor_mode %||% (if (isTRUE(cfg$use_random_censor)) "uniform" else "admin_only")
  dropout_km <- NULL
  dropout_scale_used <- 1
  if (censor_mode == "km_dropout") {
    add("\n--- 删失重构：合并两组反向KM删失分布 + 强度网格校准 ---")
    dropout_km <- pooled_dropout_km(ipd$os_trt$IPD, ipd$os_ctrl$IPD)
    add(sprintf("反向KM（删失当事件）合并中位潜在删失时间=%.2f 月（n=%d）",
                dropout_km$median, dropout_km$n))
    crn0 <- make_crn(cfg, seed, dropout_km = dropout_km, dropout_scale = 1)
    med_fu_target <- target_median_fu %||% dropout_km$median
    hr_os_target0 <- ipd_unadj_hr(ipd$os_trt$IPD, ipd$os_ctrl$IPD)
    cal_scale <- calibrate_dropout_scale(
      cfg, os_trt, os_ctrl,
      dropout_km, crn0, km_os_t, km_os_c, med_fu_target, dropout_scale_grid,
      hr_target = hr_os_target0
    )
    add(sprintf("HR校准目标：IPD未调整OS HR=%.4f（文献%.2f）", hr_os_target0, cfg$os$HR))
    add(sprintf("%-8s %8s %8s %8s %10s %10s %8s %8s", "scale", "evt_t", "evt_c", "med_FU",
                "maxΔOS_t", "maxΔOS_c", "HR", "loss"))
    for (rr in cal_scale$all) {
      add(sprintf("%-8.2f %8d %8d %8.2f %10.4f %10.4f %8.4f %8.4f",
                  rr$scale, rr$n_evt_t, rr$n_evt_c, rr$med_fu, rr$curve_t, rr$curve_c,
                  rr$hr, rr$loss))
    }
    dropout_scale_used <- cal_scale$best$scale
    add(sprintf(
      "选定 dropout_scale=%.2f（网格校准，最小化 事件比例/中位随访/OS曲线 联合误差）",
      dropout_scale_used
    ))
    add(sprintf(
      "校准后：事件 trt/ctrl=%d/%d（目标%d/%d，IPD）| 中位随访=%.2f（目标%.2f）| OS_HR=%.4f",
      cal_scale$best$n_evt_t, cal_scale$best$n_evt_c,
      cfg$os$events_trt, cfg$os$events_ctrl, cal_scale$best$med_fu, med_fu_target,
      cal_scale$best$hr
    ))
    add("说明：C_obs=min(C_admin, C_dropout)；C_dropout 来自两组合并反向KM，")
    add("      按 S_dropout(t)^scale 做强度网格校准（scale=1为原始估计强度）。")
    add(sprintf(
      "说明：本节 OS_HR 为单次CRN在试验规模(n=%d)下的实现，events≈%d/%d 时 logHR 的抽样SE≈%.3f，",
      cfg$n_trt + cfg$n_ctrl, cal_scale$best$n_evt_t, cal_scale$best$n_evt_c,
      sqrt(1 / max(cal_scale$best$n_evt_t, 1) + 1 / max(cal_scale$best$n_evt_c, 1))
    ))
    add("      与目标的偏差多在噪声范围内；治疗效应的稳健判断见下方大样本 HR_PFS_DGM 检验。\n")
  }

  # CRN + 统一 OS 验证（第1–3步指纹）
  crn <- make_crn(cfg, seed, dropout_km = dropout_km, dropout_scale = dropout_scale_used)
  D_t <- exp(os_trt$mu  + os_trt$sigma  * crn$z_D_t)
  D_c <- exp(os_ctrl$mu + os_ctrl$sigma * crn$z_D_c)
  os_met <- trial_scale_metrics(
    D_t, D_c, crn$enroll_t, crn$enroll_c, crn$C_rand_t, crn$C_rand_c,
    cfg$analysis_time
  )
  add(sprintf(
    "\n【统一验证 trial_scale_metrics】OS 事件 trt/ctrl=%d/%d (%.1f%%/%.1f%%)",
    os_met$n_os_evt_t, os_met$n_os_evt_c,
    100 * os_met$n_os_evt_t / cfg$n_trt, 100 * os_met$n_os_evt_c / cfg$n_ctrl
  ))
  add(sprintf(
    "  对照 IPD %d/%d (%.1f%%/%.1f%%)；文献 %d/%d",
    os_trt$n_events, os_ctrl$n_events,
    100 * os_trt$n_events / cfg$n_trt, 100 * os_ctrl$n_events / cfg$n_ctrl,
    cfg$os$events_trt, cfg$os$events_ctrl
  ))
  add(sprintf(
    "  回放中位 OS %.3f/%.3f | Cox OS HR=%.4f（文献 %.2f）  ← 第5步必须复现此 HR",
    os_met$med_os_t, os_met$med_os_c, os_met$HR_os, cfg$os$HR
  ))
  if (os_met$n_os_evt_t / cfg$n_trt + 0.03 < os_trt$n_events / cfg$n_trt ||
      os_met$n_os_evt_c / cfg$n_ctrl + 0.03 < os_ctrl$n_events / cfg$n_ctrl) {
    add("  警告：行政删失下事件数仍明显低于 IPD → 检查 analysis_time / enroll_max。")
  } else if (os_met$n_os_evt_t / cfg$n_trt > os_trt$n_events / cfg$n_trt + 0.05 ||
             os_met$n_os_evt_c / cfg$n_ctrl > os_ctrl$n_events / cfg$n_ctrl + 0.05) {
    add("  提示：行政删失下事件数高于 IPD → 可考虑开启随机失访并校准 cens_upper。")
  } else {
    add("  事件数与 IPD 大致接近（仅行政删失）。")
  }


  add("--- 第4步：锁定 death_par（改 R/ρ 不再修改）---")
  add(sprintf("death_par_trt  <- c(mu=%s, sigma=%.4f)", fmt_mu(os_trt$mu), os_trt$sigma))
  add(sprintf("death_par_ctrl <- c(mu=%s, sigma=%.4f)", fmt_mu(os_ctrl$mu), os_ctrl$sigma))
  if (os_trt$best_effort || os_ctrl$best_effort ||
      !os_trt$eval$pass || !os_ctrl$eval$pass) {
    add(sprintf("警告：至少一臂 OS 未达 max|ΔS|≤%.2f；已锁定最优。勿用 R/ρ 修 OS。\n", os_max_abs_tol))
  } else {
    add("OS 拟合合格，death_par 已固定。\n")
  }

  # ---- 第5步：各 ρ（本试验权重 / 偏好带 / 容差）----
  cal <- cfg$calibrate
  hr_pfs_ipd <- ipd_unadj_hr(ipd$pfs_trt$IPD, ipd$pfs_ctrl$IPD)
  hr_os_ipd  <- ipd_unadj_hr(ipd$os_trt$IPD,  ipd$os_ctrl$IPD)
  hr_crn <- make_hr_crn(cfg, seed, dropout_km = dropout_km, dropout_scale = dropout_scale_used)

  add(sprintf("--- 第5步：多目标校准（中位%s + KM + max|ΔS|(统一口径) + 容差HR）---",
              ifelse(median_mode == "soft",
                     sprintf("软约束≤%.0f%%", 100 * cal$pfs_med_rel_tol), "硬约束")))
  add("L = L_KM + w_max*(max|ΔS|_t+max|ΔS|_c) + w_H*(log HR - log target)^2 [+软中位罚]")
  add("目标函数与最终诊断的 max|ΔS| 均调用同一 evaluate_pfs_curve()、同一时间网格（无t_cap截断，")
  add(sprintf("按密度保序抽稀至≤80点），并输出有效支持区间（风险集占比≥%.0f%%）内的 max|ΔS|。",
              100 * support_frac_threshold))
  add(sprintf(
    "本试验权重 w_S=%g w_max=%g w_H=%g w_m=%g | HR目标(IPD)=%.4f | 偏好带[%.2f,%.2f] | N_HR=%d",
    cal$w_s, cal$w_max, cal$w_h, cal$w_m %||% 0, hr_pfs_ipd, cal$hr_pref_lo, cal$hr_pref_hi, cal$n_hr
  ))
  add(sprintf(
    "本试验容差：|relHR|≤%.0f%% | |rel中位|≤%.0f%% | max|ΔS|(支持区间)≤%.2f | 里程碑|ΔS|≤%.2f",
    100 * cal$hr_rel_tol, 100 * cal$pfs_med_rel_tol,
    cal$pfs_max_abs_tol, cal$pfs_milestone_tol
  ))
  add(sprintf("文献 PFS HR=%.2f（外部验证）| IPD未调整 OS HR=%.4f（文献 %.2f）",
              cfg$pfs$HR, hr_os_ipd, cfg$os$HR))
  if (!is.null(rho_main)) {
    add(sprintf("配置的主分析 ρ_main=%.2f；其他配置 ρ=%s",
                rho_main, paste(sprintf("%.2f", setdiff(rho_grid, rho_main)), collapse = ", ")))
  }
  add(sprintf("%-5s %8s %8s %8s %8s %8s %8s %8s %8s %8s",
              "rho", "mPFS_t", "mPFS_c", "HR_dgm", "HR_trl", "HR_OS",
              "relHR%", "maxΔt", "maxΔc", "L_R"))

  rho_results <- list()
  for (rho in rho_grid) {
    cat(sprintf("  rho=%.2f ...\n", rho))
    res <- calibrate_R_given_rho(
      rho, os_trt, os_ctrl,
      ipd$pfs_trt$IPD, ipd$pfs_ctrl$IPD,
      cfg, crn, hr_crn,
      med_pfs_anchor_t, med_pfs_anchor_c,
      hr_pfs_target = hr_pfs_ipd,
      pfs_milestone_months = pfs_milestone_months,
      median_mode = median_mode,
      support_frac_threshold = support_frac_threshold
    )
    rho_results[[as.character(rho)]] <- res
    add(sprintf(
      "%-5.2f %8.3f %8.3f %8.3f %8.3f %8.3f %8.1f %8.3f %8.3f %8.4f",
      rho, res$med_pfs_trt, res$med_pfs_ctrl,
      res$HR_pfs_dgm, res$HR_pfs_trial, res$HR_os,
      100 * res$hr_pfs_rel_err,
      res$pfs_diag_trt$max_abs, res$pfs_diag_ctrl$max_abs,
      res$L_R
    ))
  }

  # CRN / 统一验证一致性检查
  Dbars_t <- vapply(rho_results, `[[`, numeric(1), "D_trt_mean")
  Dbars_c <- vapply(rho_results, `[[`, numeric(1), "D_ctrl_mean")
  hr_os_all <- vapply(rho_results, `[[`, numeric(1), "HR_os")
  evt_t_all <- vapply(rho_results, `[[`, numeric(1), "n_os_evt_trt")
  add(sprintf(
    "\nCRN：潜在 D 均值跨 ρ 极差 trt=%.3e  ctrl=%.3e",
    diff(range(Dbars_t)), diff(range(Dbars_c))
  ))
  add(sprintf(
    "统一验证：各 ρ 的 OS HR = %s（应全等于第1–3步 %.4f）",
    paste(sprintf("%.4f", hr_os_all), collapse = ", "), os_met$HR_os
  ))
  add(sprintf(
    "统一验证：各 ρ 的 OS 事件(trt) = %s（应全等于 %d）",
    paste(evt_t_all, collapse = ", "), os_met$n_os_evt_t
  ))
  if (max(abs(hr_os_all - os_met$HR_os)) > 1e-10 ||
      any(evt_t_all != os_met$n_os_evt_t)) {
    add("严重错误：第5步 OS HR/事件数与第1–3步不一致 → 不可选参数。")
  } else {
    add("通过：第1–3步与第5步 OS 事件数/HR 完全一致。")
  }

  med_pfs_t_all <- vapply(rho_results, `[[`, numeric(1), "med_pfs_trt")
  med_pfs_c_all <- vapply(rho_results, `[[`, numeric(1), "med_pfs_ctrl")
  hr_pfs_dgm_all <- vapply(rho_results, `[[`, numeric(1), "HR_pfs_dgm")
  add(sprintf(
    "\n配置 ρ 间的拟合差异：PFS中位(trt)极差=%.3f  (ctrl)极差=%.3f  HR_PFS_DGM极差=%.3f",
    diff(range(med_pfs_t_all)), diff(range(med_pfs_c_all)), diff(range(hr_pfs_dgm_all))
  ))
  add("本入口按配置的 ρ 复现基线拟合，并报告拟合诊断及最终 μ,Σ。")
  add("主分析参数以锁定配置为准；ρ 无法仅由独立 PFS/OS 边际识别。")

  ready_debug <- max(abs(hr_os_all - os_met$HR_os)) <= 1e-10

  by_rho <- list()
  paper_ready_rhos <- character()

  for (rho in rho_grid) {
    r <- rho_results[[as.character(rho)]]

    tol_hr <- r$hr_pfs_rel_err <= cal$hr_rel_tol
    tol_med <- r$med_rel_err_t <= cal$pfs_med_rel_tol &&
      r$med_rel_err_c <= cal$pfs_med_rel_tol
    tol_max <- r$pfs_diag_trt$max_abs_supported <= cal$pfs_max_abs_tol &&
      r$pfs_diag_ctrl$max_abs_supported <= cal$pfs_max_abs_tol
    tol_ms <- r$pfs_diag_trt$max_milestone <= cal$pfs_milestone_tol &&
      r$pfs_diag_ctrl$max_milestone <= cal$pfs_milestone_tol
    ready_paper <- ready_debug && tol_hr && tol_med && tol_max && tol_ms
    if (ready_paper) paper_ready_rhos <- c(paper_ready_rhos, sprintf("%.2f", rho))

    reasons <- character()
    if (!ready_debug) reasons <- c(reasons, "OS HR 第1–3步与第5步不一致")
    if (!tol_hr) {
      reasons <- c(reasons, sprintf("HR rel误差%.1f%% > %.0f%% (HR=%.3f)",
                                   100 * r$hr_pfs_rel_err, 100 * cal$hr_rel_tol, r$HR_pfs_dgm))
    }
    if (!tol_med) {
      reasons <- c(reasons, sprintf("PFS中位rel误差 %.1f%%/%.1f%% > %.0f%%",
                                   100 * r$med_rel_err_t, 100 * r$med_rel_err_c,
                                   100 * cal$pfs_med_rel_tol))
    }
    if (!tol_max) {
      reasons <- c(reasons, sprintf("支持区间max|ΔS|=%.3f/%.3f > %.2f",
                                   r$pfs_diag_trt$max_abs_supported,
                                   r$pfs_diag_ctrl$max_abs_supported,
                                   cal$pfs_max_abs_tol))
    }
    if (!tol_ms) {
      reasons <- c(reasons, sprintf("里程碑max|ΔS|=%.3f/%.3f > %.2f",
                                   r$pfs_diag_trt$max_milestone,
                                   r$pfs_diag_ctrl$max_milestone,
                                   cal$pfs_milestone_tol))
    }

    Sig_t <- make_Sigma(r$sig_R_trt, os_trt$sigma, rho)
    Sig_c <- make_Sigma(r$sig_R_ctrl, os_ctrl$sigma, rho)
    mu_t <- c(r$mu_R_trt, os_trt$mu)
    mu_c <- c(r$mu_R_ctrl, os_ctrl$mu)

    add(sprintf("\n===== ρ=%.2f%s =====", rho,
                ifelse(!is.null(rho_main) && isTRUE(all.equal(rho, rho_main)), "（主分析）", "（配置ρ拟合）")))
    add(sprintf("  μ_R,σ_R trt  : %s , %.4f  | 理论 S_PFS(m)=%.4f",
                fmt_mu(r$mu_R_trt), r$sig_R_trt, r$S_pfs_at_m_trt))
    add(sprintf("  μ_R,σ_R ctrl : %s , %.4f  | 理论 S_PFS(m)=%.4f",
                fmt_mu(r$mu_R_ctrl), r$sig_R_ctrl, r$S_pfs_at_m_ctrl))
    add(sprintf(
      "  大样本 DGM HR_PFS=%.4f | 目标=%.4f | rel误差=%.1f%% %s | 文献=%.2f | 偏好带[%g,%g]%s",
      r$HR_pfs_dgm, hr_pfs_ipd, 100 * r$hr_pfs_rel_err,
      ifelse(tol_hr, "PASS", "FAIL"),
      cfg$pfs$HR, cal$hr_pref_lo, cal$hr_pref_hi,
      ifelse(r$hr_in_pref_band, " ✓", "")
    ))
    add(sprintf(
      "  试验规模回放 HR_PFS=%.4f  HR_OS=%.4f（应=第1–3步 %.4f）",
      r$HR_pfs_trial, r$HR_os, os_met$HR_os
    ))
    add(sprintf(
      "  回放中位 PFS %.3f/%.3f（锚定 %.3f/%.3f；rel误差 %.1f%%/%.1f%%）%s",
      r$med_pfs_trt, r$med_pfs_ctrl,
      med_pfs_anchor_t, med_pfs_anchor_c,
      100 * r$med_rel_err_t, 100 * r$med_rel_err_c,
      ifelse(tol_med, " PASS", " FAIL")
    ))
    add(sprintf("  OS 中位 %.3f/%.3f | OS 事件 %d/%d",
                r$med_os_trt, r$med_os_ctrl,
                r$n_os_evt_trt, r$n_os_evt_ctrl))
    add(sprintf("  L_R=%.4f (D_KM=%.4f+%.4f, max|ΔS|_opt=%.3f+%.3f, HR罚=%.4f)",
                r$L_R, r$D_km_t, r$D_km_c,
                r$max_abs_opt_t, r$max_abs_opt_c, r$hr_pen))

    add("  PFS 曲线诊断（理论 DGM vs IPD KM；目标函数与诊断统一口径）：")
    for (lab in c("trt", "ctrl")) {
      dg <- if (lab == "trt") r$pfs_diag_trt else r$pfs_diag_ctrl
      add(sprintf(
        "    [%s] S(m)=%.4f | max_t|ΔS|(全程)=%.4f  @t=%.2fm (KM=%.4f DGM=%.4f, 风险集=%d/%.0f%%) | 里程碑max|ΔS|=%.4f %s",
        lab, dg$S_at_m, dg$max_abs, dg$t_at_max, dg$S_km_at_max, dg$S_model_at_max,
        dg$n_risk_at_max, 100 * dg$frac_risk_at_max,
        dg$max_milestone,
        ifelse(dg$max_milestone <= cal$pfs_milestone_tol, "PASS", "FAIL")
      ))
      add(sprintf(
        "        支持区间(风险集≥%.0f%%)内 max|ΔS|=%.4f %s  @t=%.2fm (KM=%.4f DGM=%.4f, 风险集=%d/%.0f%%)%s",
        100 * support_frac_threshold, dg$max_abs_supported,
        ifelse(dg$max_abs_supported <= cal$pfs_max_abs_tol, "PASS", "FAIL"),
        dg$t_at_max_supported, dg$S_km_at_max_supported, dg$S_model_at_max_supported,
        dg$n_risk_at_max_supported, 100 * dg$frac_risk_at_max_supported,
        ifelse(dg$max_in_thin_tail,
               "  [全程最大误差发生在风险集稀薄尾部→按支持区间口径评估]",
               "  [全程最大误差本身即在风险集充足区间→提示单一Lognormal可能失配]")
      ))
      for (i in seq_len(nrow(dg$milestones))) {
        m <- dg$milestones[i, ]
        add(sprintf(
          "         %2dm: 模型=%.3f  KM=%.3f  |Δ|=%.3f %s",
          m$month, m$S_model, m$S_km, m$abs_err,
          ifelse(m$abs_err <= cal$pfs_milestone_tol, "PASS", "FAIL")
        ))
      }
    }

    add(sprintf("  mu.trt  <- %s", fmt_pair(mu_t)))
    add(sprintf("  mu.ctrl <- %s", fmt_pair(mu_c)))
    add(sprintf(
      "  Sigma.trt  <- matrix(c(%.4f, %.4f, %.4f, %.4f), nrow=2, byrow=TRUE)  # rho=%.2f",
      Sig_t[1, 1], Sig_t[1, 2], Sig_t[2, 1], Sig_t[2, 2], rho
    ))
    add(sprintf(
      "  Sigma.ctrl <- matrix(c(%.4f, %.4f, %.4f, %.4f), nrow=2, byrow=TRUE)",
      Sig_c[1, 1], Sig_c[1, 2], Sig_c[2, 1], Sig_c[2, 2]
    ))
    add("  【多目标容差检查】")
    add(sprintf("    |rel HR|≤%.0f%%          : %s (%.1f%%, HR=%.3f)",
                100 * cal$hr_rel_tol, ifelse(tol_hr, "PASS", "FAIL"),
                100 * r$hr_pfs_rel_err, r$HR_pfs_dgm))
    add(sprintf("    |rel PFS中位|≤%.0f%%     : %s (%.1f%%/%.1f%%)",
                100 * cal$pfs_med_rel_tol, ifelse(tol_med, "PASS", "FAIL"),
                100 * r$med_rel_err_t, 100 * r$med_rel_err_c))
    add(sprintf("    支持区间max_t|ΔS|≤%.2f  : %s (%.3f/%.3f；全程 %.3f/%.3f)",
                cal$pfs_max_abs_tol, ifelse(tol_max, "PASS", "FAIL"),
                r$pfs_diag_trt$max_abs_supported, r$pfs_diag_ctrl$max_abs_supported,
                r$pfs_diag_trt$max_abs, r$pfs_diag_ctrl$max_abs))
    add(sprintf("    里程碑|ΔS|≤%.2f       : %s (%.3f/%.3f)",
                cal$pfs_milestone_tol, ifelse(tol_ms, "PASS", "FAIL"),
                r$pfs_diag_trt$max_milestone, r$pfs_diag_ctrl$max_milestone))
    add(sprintf("  【全部拟合检查】%s", ifelse(ready_paper, "PASS", "FAIL")))
    if (length(reasons)) add(paste0("    未达标：", paste(reasons, collapse = "；")))

    by_rho[[sprintf("%.2f", rho)]] <- list(
      rho = rho,
      mu.trt = mu_t, mu.ctrl = mu_c,
      Sigma.trt = Sig_t, Sigma.ctrl = Sig_c,
      HR_pfs_dgm = r$HR_pfs_dgm, HR_pfs_trial = r$HR_pfs_trial,
      med_pfs_trt = r$med_pfs_trt, med_pfs_ctrl = r$med_pfs_ctrl,
      hr_pfs_rel_err = r$hr_pfs_rel_err,
      med_rel_err_t = r$med_rel_err_t, med_rel_err_c = r$med_rel_err_c,
      max_abs_t = r$pfs_diag_trt$max_abs, max_abs_c = r$pfs_diag_ctrl$max_abs,
      max_abs_supported_t = r$pfs_diag_trt$max_abs_supported,
      max_abs_supported_c = r$pfs_diag_ctrl$max_abs_supported,
      max_milestone_t = r$pfs_diag_trt$max_milestone,
      max_milestone_c = r$pfs_diag_ctrl$max_milestone,
      tol_hr = tol_hr, tol_med = tol_med, tol_max = tol_max, tol_ms = tol_ms,
      ready_for_paper_main_sim = ready_paper,
      is_main = !is.null(rho_main) && isTRUE(all.equal(rho, rho_main))
    )
  }

  add(sprintf(
    "\n满足全部指定拟合容差的 ρ：%s",
    ifelse(length(paper_ready_rhos), paste(paper_ready_rhos, collapse = ", "), "无")
  ))
  add(sprintf("\n【OS HR一致性检查】%s", ifelse(ready_debug, "PASS", "FAIL")))
  add(sprintf("说明：本试验（%s）入组/删失/权重/偏好带/容差均在本试验专属脚本中独立配置，", trial_name))
  add("      各试验的拟合设置独立；当前入口报告配置 ρ 的基线拟合，主分析采用锁定配置。")

  writeLines(report, file.path(out_dir, "calibration_report.txt"))

  sink(file.path(out_dir, "params_for_data_simulation.R"))
  cat(sprintf("# OS-first | anchor=%s | %s | rho_grid=%s (按配置的ρ输出基线拟合结果)\n",
              anchor_source, trial_name, paste(sprintf("%.2f", rho_grid), collapse = ",")))
  cat(sprintf("# PFS HR target = IPD unadj %.4f (lit %.2f external only)\n",
              hr_pfs_ipd, cfg$pfs$HR))
  cat(sprintf("# censor: analysis_time=%g, enroll_max=%g, random=%s",
              cfg$analysis_time, cfg$enroll_max, isTRUE(cfg$use_random_censor)))
  if (isTRUE(cfg$use_random_censor)) cat(sprintf(", cens_upper=%g", cfg$cens_upper))
  cat("\n")
  cat(sprintf("# by_rho[[\"%.2f\"]]$mu.trt 等：按配置的 ρ 读取基线参数；主分析参数以锁定配置为准\n", rho_grid[1]))
  cat(sprintf("%s_os_first <- list(\n", gsub("[^A-Za-z0-9_]", "_", trial_name)))
  cat(sprintf("  rho_grid = c(%s),\n", paste(sprintf("%.2f", rho_grid), collapse = ", ")))
  cat(sprintf("  paper_ready_rhos = c(%s),\n",
              ifelse(length(paper_ready_rhos),
                     paste(sprintf("\"%s\"", paper_ready_rhos), collapse = ", "), "")))
  cat("  by_rho = list(\n")
  for (rho in rho_grid) {
    key <- sprintf("%.2f", rho)
    b <- by_rho[[key]]
    cat(sprintf("    `%s` = list(\n", key))
    cat(sprintf("      rho = %.2f,\n", b$rho))
    cat(sprintf("      mu.trt = %s,\n", fmt_pair(b$mu.trt)))
    cat(sprintf("      mu.ctrl = %s,\n", fmt_pair(b$mu.ctrl)))
    cat(sprintf("      Sigma.trt = matrix(c(%.4f, %.4f, %.4f, %.4f), nrow=2, byrow=TRUE),\n",
                b$Sigma.trt[1, 1], b$Sigma.trt[1, 2], b$Sigma.trt[2, 1], b$Sigma.trt[2, 2]))
    cat(sprintf("      Sigma.ctrl = matrix(c(%.4f, %.4f, %.4f, %.4f), nrow=2, byrow=TRUE),\n",
                b$Sigma.ctrl[1, 1], b$Sigma.ctrl[1, 2], b$Sigma.ctrl[2, 1], b$Sigma.ctrl[2, 2]))
    cat(sprintf("      HR_pfs_dgm = %.6f,\n", b$HR_pfs_dgm))
    cat(sprintf("      med_pfs = c(%.4f, %.4f),\n", b$med_pfs_trt, b$med_pfs_ctrl))
    cat(sprintf("      max_abs_S = c(%.4f, %.4f),\n", b$max_abs_t, b$max_abs_c))
    cat(sprintf("      max_abs_S_supported = c(%.4f, %.4f),\n",
                b$max_abs_supported_t, b$max_abs_supported_c))
    cat(sprintf("      is_main = %s,\n", b$is_main))
    cat(sprintf("      ready_for_paper_main_sim = %s\n", b$ready_for_paper_main_sim))
    cat(sprintf("    )%s\n", ifelse(rho == rho_grid[length(rho_grid)], "", ",")))
  }
  cat("  ),\n")
  cat(sprintf("  analysis_time = %g,\n", cfg$analysis_time))
  cat(sprintf("  enroll_max = %g,\n", cfg$enroll_max))
  cat(sprintf("  use_random_censor = %s,\n", isTRUE(cfg$use_random_censor)))
  cat(sprintf("  cens_upper = %g,\n", cfg$cens_upper))
  cat(sprintf("  hr_pref = c(%.4f, %.4f),\n", cal$hr_pref_lo, cal$hr_pref_hi))
  cat(sprintf("  death_par_trt = c(mu=%.6f, sigma=%.6f),\n", os_trt$mu, os_trt$sigma))
  cat(sprintf("  death_par_ctrl = c(mu=%.6f, sigma=%.6f),\n", os_ctrl$mu, os_ctrl$sigma))
  cat(sprintf("  HR_pfs_target_ipd = %.6f,\n", hr_pfs_ipd))
  cat(sprintf("  HR_pfs_lit = %.4f,\n", cfg$pfs$HR))
  cat(sprintf("  os_HR_trial = %.6f,\n", os_met$HR_os))
  cat(sprintf("  os_events_trial = c(%d, %d),\n", os_met$n_os_evt_t, os_met$n_os_evt_c))
  cat(sprintf("  ready_for_debug_sim = %s\n)\n", ready_debug))
  sink()

  saveRDS(list(
    anchor = anchor_source,
    trial = trial_name,
    cfg_censor = list(
      analysis_time = cfg$analysis_time,
      enroll_max = cfg$enroll_max,
      use_random_censor = isTRUE(cfg$use_random_censor),
      cens_upper = cfg$cens_upper
    ),
    cfg_calibrate = cal,
    os_trt = os_trt, os_ctrl = os_ctrl,
    os_met_step13 = os_met,
    hr_pfs_ipd = hr_pfs_ipd, hr_os_ipd = hr_os_ipd,
    med_os_anchor = c(trt = med_os_anchor_t, ctrl = med_os_anchor_c),
    med_pfs_anchor = c(trt = med_pfs_anchor_t, ctrl = med_pfs_anchor_c),
    rho_results = rho_results, by_rho = by_rho,
    paper_ready_rhos = paper_ready_rhos,
    ready_for_debug_sim = ready_debug
  ), file.path(out_dir, "calibrated_params.rds"))

  cat(paste(report, collapse = "\n"))
  cat(sprintf("\n→ %s\n", file.path(out_dir, "calibration_report.txt")))
  cat(sprintf("\n完成：%s。结果在 %s/\n", trial_name, out_dir))

  invisible(list(report = report, out_dir = out_dir))
}

# ------------------------------------------------------------

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
cfg <- list(
  trial_dir = .script_dir,
  output_dir = .output_dir,
  rds = file.path(.package_dir, "inputs/reconstructed_data/fitting_anchors/impower131/IMpower131_marginal_ipd.rds"),
  n_trt = 343, n_ctrl = 340,
  # --- 入组 / 删失（IMpower 专属）---
  analysis_time = 41,
  enroll_max = 22,
  use_random_censor = FALSE,
  cens_upper = 300,          # 仅 use_random_censor=TRUE 时启用
  # --- 文献报告（外部验证）---
  pfs = list(median_trt = 6.3, median_ctrl = 5.6, HR = 0.71,
             events_trt = 270, events_ctrl = 289),
  os  = list(median_trt = 14.2, median_ctrl = 13.5, HR = 0.88,
             events_trt = 228, events_ctrl = 245),
  # --- 本试验校准超参（IMpower 专属，不影响 OAK）---
  calibrate = list(
    w_s = 1, w_max = 4, w_h = 12, n_hr = 8000L,
    hr_pref_lo = 0.70, hr_pref_hi = 0.74,   # 围绕文献/IPD ~0.70–0.72
    hr_rel_tol = 0.05,
    pfs_med_rel_tol = 0.05,
    pfs_max_abs_tol = 0.05,
    pfs_milestone_tol = 0.05
  )
)

run_os_first_calibration(
  trial_name = "IMpower131",
  cfg = cfg,
  rho_grid = 0.5,  # adopted E05-v6 baseline only
  anchor_source = "ipd",
  out_subdir = "os_first_method_ipd_anchor",
  seed = 2026L,
  os_max_abs_tol = 0.03,
  os_milestone_months = c(6, 12, 18),
  pfs_milestone_months = c(6, 12, 18),
  mu_d_med_tol = 0.02
)
