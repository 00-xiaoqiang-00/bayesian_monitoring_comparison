#' OCC.Table.LG
#'
#' Bayesian monitoring with a dynamic interim futility boundary that multiplies
#' the final threshold c_F into the interim shape (n/N.max)^gamma_E:
#'
#'   Futility(n)    = c_F * (n / N.max)^gamma_E
#'   Superiority(n) = 1 - (n / N.max)^gamma_E   [optional, see enable_superiority]
#'   Final cutoff   = c_F
#'
#' At every interim look (n < N.max): stop for futility if
#' P(WR > eta_interim | data) <= Futility(n); stop early for superiority (go)
#' if PP_n >= Superiority(n); otherwise continue.
#' At the final look (n = N.max): declare "go" if
#' P(WR > eta_final | data) > c_F, else "no_go".
#' There is no continuation option at the final look.
#'
#' This is a drop-in sibling of OCC.Table(): the posterior-probability
#' computation for each design (1, 2/3/4, 6, 7) is identical to
#' OCandCutoff_final7.R; only the stopping-boundary / decision logic differs.
#'
#' @useDynLib BayesianWinRatio, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
#' @param myData A matrix with N.max rows and 7 columns (see OCC.Table).
#' @param N.max Maximum number of patients to enroll.
#' @param design 1, 2, 3, 4, 6, or 7 (same meaning as OCC.Table).
#' @param cohort Numeric vector of interim look sample sizes. N.max is
#'   automatically appended as the final look if not already the last element.
#' @param recruit.int Recruitment interval.
#' @param m0,L0,v0,S0 Prior parameters (design-specific, see OCC.Table).
#' @param time_max Upper truncation time for latent sampling.
#' @param eta Default effect-size margin for both interim and final looks if
#'   \code{eta_interim} / \code{eta_final} are not supplied.
#' @param eta_interim Optional margin for interim monitoring
#'   (e.g. \code{1} tests WR > 1 for futility screening).
#' @param eta_final Optional margin for the final go/no-go test
#'   (e.g. \code{1} or \code{1.5}).
#' @param final_analysis_time Optional common-data-cutoff time measured from
#'   first enrollment. If supplied, the final look uses this calendar time;
#'   interim looks remain enrollment-driven. If \code{NULL}, the final look
#'   occurs when enrollment reaches \code{N.max}, preserving compatibility behavior.
#' @param gamma_E Interim boundary shape parameter in
#'   Futility(n) = c_F * (n / N.max)^gamma_E.
#' @param c_F Shared strength parameter: final go requires
#'   P(WR > eta_final | data) > c_F, and the same c_F scales the interim
#'   futility boundary.
#' @param thin_MCMC Thinning.
#' @param Niter Number of MCMC iterations.
#' @param Sigma.ctrl Control variance/covariance input (design-specific).
#' @param enable_superiority Logical, default FALSE. If FALSE, the interim
#'   superiority (early efficacy) boundary is disabled (sup_bound = Inf), so
#'   the design only monitors for futility at interim looks; the trial can
#'   still declare "go" at the final analysis (n = N.max) via c_F.
#'
#' @returns A list with:
#' \tabular{ll}{
#'   \code{stop_look} \tab Index into \code{cohort} at which the trial stopped. \cr
#'   \code{n_final} \tab Sample size at the look where the trial stopped. \cr
#'   \code{is_final_look} \tab TRUE if the trial ran to N.max. \cr
#'   \code{decision} \tab One of "futility", "superiority", "go", "no_go". \cr
#'   \code{go} \tab Logical; TRUE for "superiority" or "go", FALSE otherwise. \cr
#'   \code{pts.stop} \tab Realized sample size (== n_final). \cr
#'   \code{probs} \tab Posterior probability PP at the stopping look. \cr
#'   \code{fut_bound}, \code{sup_bound} \tab Interim boundary values at the
#'     stopping look. \cr
#'   \code{c_F}, \code{final_bound} \tab Final go threshold (equal at the
#'     final look). \cr
#'   \code{eta_interim}, \code{eta_final} \tab Margins used at interim/final. \cr
#'   \code{WR} \tab Posterior mean win ratio (design 1 only; NA otherwise). \cr
#' }
#' @examples
#' \dontrun{
#' OCC.Table.LG(myData = data[,,1], N.max = 80, design = 1,
#'   cohort = c(32, 48, 64), recruit.int = 0.25,
#'   m0 = c(0,0), L0 = diag(10^6, 2), v0 = 4, S0 = diag(10^(-6), 2),
#'   time_max = 100, eta = 1, gamma_E = 0.90, c_F = 0.90,
#'   thin_MCMC = 5, Niter = 30000,
#'   Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = TRUE))
#' }

#' @export
OCC.Table.LG <- function(myData, N.max, design, cohort, recruit.int,
                         m0, L0, v0, S0, time_max, eta, gamma_E, c_F,
                         thin_MCMC, Niter, Sigma.ctrl, enable_superiority = FALSE,
                         eta_interim = NULL, eta_final = NULL,
                         final_analysis_time = NULL) {
  WR <- NA
  probs <- NA_real_
  fut_bound <- NA_real_
  sup_bound <- NA_real_
  eta_interim <- if (is.null(eta_interim)) eta else eta_interim
  eta_final <- if (is.null(eta_final)) eta else eta_final

  tryCatch({
    cohort <- sort(unique(as.integer(cohort)))
    cohort <- cohort[cohort > 0 & cohort <= N.max]
    if (length(cohort) == 0 || cohort[length(cohort)] != N.max) {
      cohort <- c(cohort, N.max)
    }

    Time.entry <- Time.current <- n.current <- j.cohort <- 0
    trial.stop <- 0
    decision <- NA_character_
    go <- NA

    while (trial.stop != 1 && j.cohort < length(cohort)) {
      j.cohort <- j.cohort + 1
      n.current <- cohort[j.cohort]
      is_final_look <- (n.current == N.max)

      Time.entry <- recruit.int * (c(1:n.current) - 1)
      enrollment_time <- recruit.int * (n.current - 1)
      Time.current <- if (is_final_look && !is.null(final_analysis_time)) {
        final_analysis_time
      } else {
        enrollment_time
      }
      if (Time.current < enrollment_time) {
        stop("Analysis time precedes enrollment of the current cohort: ",
             "Time.current=", Time.current,
             ", required>=", enrollment_time)
      }
      temp <- myData[1:n.current - 1, ]

      eta_use <- if (is_final_look) eta_final else eta_interim

      if (design == 1) {
        currentData <- temp
        result <- winratio(
          currentdd = currentData,
          n_current = n.current,
          Time_current = Time.current,
          Time_entry = Time.entry,
          m0 = m0,
          L0 = L0,
          v0 = v0,
          S0 = S0,
          time_max = time_max,
          eta = eta_use,
          lambda = 0,          # unused inside winratio() other than for its own $cutoff, which we discard
          N.max = N.max,
          thin_MCMC = thin_MCMC,
          Niter = Niter,
          Sigma.ctrl = Sigma.ctrl
        )
        probs <- result$probs
        WR <- result$WR
      } else if (design %in% c(2, 3, 4)) {
        temp[, 3] <- pmin(Time.current - Time.entry[-n.current], temp[, 3], temp[, 2])
        if (design %in% c(2, 3)) {
          currentData <- temp[, c(design - 1, 3:6)]
        } else {
          currentData <- data.frame(evnt.time = pmin(temp[, 1], temp[, 2]), temp[, c(3:6)])
          currentData <- as.matrix(currentData)
        }

        currentData[, 5] <- as.integer(currentData[, 1] <= currentData[, 2])
        currentData[, 1] <- ifelse(currentData[, 5] == 1, currentData[, 1], NA)

        n.current.trt <- sum(currentData[, 3])
        n.current.ctrl <- n.current - 1 - n.current.trt
        currentData.trt <- currentData[currentData[, 3] == 1, ]
        currentData.ctrl <- currentData[currentData[, 3] == 0, ]

        temp_var <- stats::var(log(currentData.trt[, 1]), na.rm = TRUE)
        init_trt_Sigma <- if (is.na(temp_var)) Sigma.ctrl else temp_var
        trt.post <- update_theta_univariate(
          N_iter = Niter, dd = currentData.trt, n = n.current.trt,
          L0 = L0, m0 = m0, v0 = v0, S0 = S0, time_max = time_max,
          init_Sigma = init_trt_Sigma
        )

        temp_var <- stats::var(log(currentData.ctrl[, 1]), na.rm = TRUE)
        init_ctrl_Sigma <- if (is.na(temp_var)) Sigma.ctrl else temp_var
        ctrl.post <- update_theta_univariate(
          N_iter = Niter, dd = currentData.ctrl, n = n.current.ctrl,
          L0 = L0, m0 = m0, v0 = v0, S0 = S0, time_max = time_max,
          init_Sigma = init_ctrl_Sigma
        )

        burn_MCMC <- as.integer(0.3 * Niter)
        idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
        probs <- mean(trt.post$theta[idxs, 1] - ctrl.post$theta[idxs, 1] > eta_use)
      } else if (design == 6) {
        temp[, 3] <- pmin(Time.current - Time.entry[-n.current], temp[, 3], temp[, 2])
        event_time <- pmin(temp[, 1], temp[, 2])
        censor_time <- temp[, 3]
        obs_time <- pmin(event_time, censor_time)
        status <- as.integer(event_time <= censor_time)
        group <- as.integer(temp[, 4])

        beta_samples <- .bayes_cox_mcmc_beta(
          time = obs_time, status = status, group = group,
          Niter = Niter, thin_MCMC = thin_MCMC,
          prior_mean = 0, prior_sd = 1.5, proposal_sd = 0.12
        )

        probs <- if (all(is.na(beta_samples))) NA_real_ else mean(beta_samples < 0, na.rm = TRUE)
        WR <- NA
      } else if (design == 7) {
        temp[, 3] <- pmin(Time.current - Time.entry[-n.current], temp[, 3], temp[, 2])

        currentData <- data.frame(evnt.time = pmin(temp[, 1], temp[, 2]), temp[, c(3:6)])
        currentData <- as.matrix(currentData)
        currentData[, 5] <- as.integer(currentData[, 1] <= currentData[, 2])
        currentData[, 1] <- ifelse(currentData[, 5] == 1, currentData[, 1], NA)

        n.current.trt <- sum(currentData[, 3])
        n.current.ctrl <- n.current - 1 - n.current.trt
        currentData.trt <- currentData[currentData[, 3] == 1, ]
        currentData.ctrl <- currentData[currentData[, 3] == 0, ]

        m0_u <- .as_scalar(m0, 0)
        L0_u <- .as_scalar(L0, 1e6)
        v0_u <- .as_scalar(v0, 4)
        S0_u <- .as_scalar(S0, 1e-6)
        Sigma_u <- .as_scalar(Sigma.ctrl, 1)

        temp_var <- stats::var(log(currentData.trt[, 1]), na.rm = TRUE)
        init_trt_Sigma <- if (is.na(temp_var)) Sigma_u else temp_var
        trt.post <- update_theta_univariate(
          N_iter = Niter, dd = currentData.trt, n = n.current.trt,
          L0 = L0_u, m0 = m0_u, v0 = v0_u, S0 = S0_u,
          time_max = time_max, init_Sigma = init_trt_Sigma
        )

        temp_var <- stats::var(log(currentData.ctrl[, 1]), na.rm = TRUE)
        init_ctrl_Sigma <- if (is.na(temp_var)) Sigma_u else temp_var
        ctrl.post <- update_theta_univariate(
          N_iter = Niter, dd = currentData.ctrl, n = n.current.ctrl,
          L0 = L0_u, m0 = m0_u, v0 = v0_u, S0 = S0_u,
          time_max = time_max, init_Sigma = init_ctrl_Sigma
        )

        burn_MCMC <- as.integer(0.3 * Niter)
        idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
        tau <- Time.current

        trt_rmst <- mapply(.rmst_lognormal, trt.post$theta[idxs, 1], trt.post$theta[idxs, 2],
                           MoreArgs = list(tau = tau))
        ctrl_rmst <- mapply(.rmst_lognormal, ctrl.post$theta[idxs, 1], ctrl.post$theta[idxs, 2],
                            MoreArgs = list(tau = tau))

        probs <- mean((trt_rmst - ctrl_rmst) > 0, na.rm = TRUE)
        WR <- NA
      } else {
        stop("Unsupported design value. Supported: 1,2,3,4,6,7")
      }

      ratio <- n.current / N.max
      fut_bound <- c_F * (ratio^gamma_E)
      sup_bound <- if (isTRUE(enable_superiority)) {
        1 - ratio^gamma_E
      } else {
        Inf # superiority stopping disabled: never triggers at interim looks
      }

      if (!is_final_look) {
        if (!is.na(probs) && probs <= fut_bound) {
          trial.stop <- 1
          decision <- "futility"
          go <- FALSE
        } else if (!is.na(probs) && probs >= sup_bound) {
          trial.stop <- 1
          decision <- "superiority"
          go <- TRUE
        }
        # otherwise: continue to next look
      } else {
        trial.stop <- 1
        if (!is.na(probs) && probs > c_F) {
          decision <- "go"
          go <- TRUE
        } else {
          decision <- "no_go"
          go <- FALSE
        }
      }
    }

    list(
      stop_look = j.cohort,
      n_final = n.current,
      is_final_look = (n.current == N.max),
      decision = decision,
      go = go,
      pts.stop = n.current,
      probs = probs,
      fut_bound = fut_bound,
      sup_bound = sup_bound,
      final_bound = if (is_final_look) c_F else NA_real_,
      gamma_E = gamma_E,
      c_F = c_F,
      eta_interim = eta_interim,
      eta_final = eta_final,
      WR = WR
    )
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
  })
}
