# Case config: IMpower131
#
# Parameters correspond to the locked IMpower131 design configuration.
# (calibrate_LG.R + validate_LG_scenarios.R + H1_8_18 / observed caches).
# Do not edit global scripts for this case; edit this file only.
#
# NOTE on mu:
#   Code uses bivariate lognormal: T = exp(X), X ~ N(mu, Sigma).
#   With mu = log(a), Median(T) = a; Mean(T) = a * exp(sigma^2/2).

case_id   <- "IMpower131"
case_note <- paste(
  "IMpower131 LG win-ratio design; locked common DGM E05 v6;",
  "near-administrative censoring plus a common 41-month calendar cutoff"
)
dgm_version <- "IMpower131_bivariate_lognormal_v6_E05_locked_20260801"

# ---- trial design skeleton ----
N.max          <- 680
interim_cohort <- c(272L, 408L, 544L)   # 40%, 60%, 80%; N.max is final
recruit.int    <- 22 / 680
cens_upper     <- 1e6
final_analysis_time <- 41                # common data cutoff, months from first enrollment

# ---- MCMC / simulation ----
design       <- 1L
m0           <- c(0, 0)
L0           <- diag(10^6, 2)
v0           <- 4
S0           <- diag(10^(-6), 2)
time_max     <- 100
thin_MCMC    <- 5
Niter        <- 10000
L_sim <- 10000
mc.cores <- if (.Platform$OS.type == "windows") 1L else 70L
seed <- 26073001L

# ---- decision / calibration targets ----
eta_interim        <- 1
eta_final          <- 1
enable_superiority <- FALSE
alpha <- 0.05
power_target       <- 0.80
# Grid used in later IMpower exploration (c_F from 0.85)
c_F_grid <- seq(0.80, 0.99, by = 0.005)
gamma_E_grid <- seq(0.00, 4.00, by = 0.05)
# FALSE: reuse each compatible cache; build only missing/incompatible ones.
# TRUE: force-rebuild every cache (escape hatch; delete .rds for a single redo).
rebuild_pp_cache   <- FALSE

# Recommended fixed design from IMpower analysis (for validation replay)
c_F_recommended     <- 0.85
gamma_E_recommended <- 1.00

# ---- covariance (log-time) ----
Sigma.ctrl <- matrix(
  c(0.4209, 0.3287,
    0.3287, 1.0268),
  nrow = 2, byrow = TRUE
)
Sigma_h0_trt <- Sigma.ctrl
Sigma_h1_trt <- Sigma.ctrl

# ---- scenario means on LOG scale: mu = E[log T] ----
# Control / H0: identical arms, exp(mu) = 6.1202, 13.4685.
# The first value is the latent progression-time median; observed PFS is min(R,D).
mu_h0_trt  <- c(log(6.1202), log(13.4685))
mu_h0_ctrl <- c(log(6.1202), log(13.4685))

# Calibration alternative: treatment exp(mu)=8,18; control exp(mu)=6.1202,13.4685
mu_h1_trt  <- c(log(8), log(18))
mu_h1_ctrl <- c(log(6.1202), log(13.4685))

# Extra validation scenario from observed treatment-arm data
scenarios_extra <- list(
  list(
    name = "H1_observed",
    role = "alternative",
    note = "E05 v6 observed treatment exp(mu)=7.5610,14.3431",
    cache_file = "calib_stats_H1_observed.rds",
    mu.trt = c(log(7.5610), log(14.3431)),
    mu.ctrl = c(log(6.1202), log(13.4685)),
    Sigma.trt = matrix(
      c(0.8207, 0.5298,
        0.5298, 1.3679),
      nrow = 2, byrow = TRUE
    ),
    Sigma.ctrl = Sigma.ctrl
  )
)

# Global output/cache/* were built with older scenario parameters.
# They are NOT compatible with this case configuration.
link_global_caches <- FALSE

case_root <- NULL

# ---- analysis schedule version ----
design_version <- "calibration_final_all_n_v3"
analysis_row_rule <- "interim_n_minus_1_final_all_n_v2"
interim_fractions <- c(0.40, 0.60, 0.80)
parent_project <- ""

# Formal calibration provenance. Decision parameters are selected only after
# all H0/H1 PP caches are complete.
calibration_role <- "formal_parameter_calibration"
calibration_outer_replicates <- 10000L
calibration_final_rule <- "PP_final > c_F"
calibration_futility_rule <- "PP_k <= b_k"
calibration_selection_constraints <- paste(
  "0.04 <= T1E <= 0.05;",
  "Power >= 0.80;",
  "H1 PET <= 0.10;",
  "c_F >= 0.82"
)
