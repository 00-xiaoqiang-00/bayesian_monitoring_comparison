# Statistical functions for exact-count screening and paired-MCSE selection.
# The command-line drivers restrict their scope to Cox and RMST.

prepare_cache <- function(path) {
  x <- readRDS(path)
  attempted <- length(x)
  if (attempted != 10000L) stop("Expected 10,000 records: ", path)
  valid <- vapply(x, function(z) {
    flag <- if (!is.null(z$ok)) isTRUE(z$ok) else if (!is.null(z$fit_ok)) isTRUE(z$fit_ok) else TRUE
    pp <- if (!is.null(z$PP)) z$PP else z$PP_interims
    flag && length(pp) == 3L && all(is.finite(pp)) &&
      length(z$PP_final) == 1L && is.finite(z$PP_final)
  }, logical(1L))
  if (mean(valid) < 0.95) stop("Valid-fit rate below 95%: ", path)
  y <- x[valid]
  list(
    pp = do.call(rbind, lapply(y, function(z) as.numeric(if (!is.null(z$PP)) z$PP else z$PP_interims))),
    final = vapply(y, function(z) as.numeric(z$PP_final), numeric(1L)),
    attempted = attempted,
    valid = sum(valid),
    failed = attempted - sum(valid)
  )
}

scenario_metrics <- function(cache, pp_threshold, lambda_prime, interim_counts, n_max) {
  boundary <- pp_threshold * (interim_counts / n_max)^lambda_prime
  below <- sweep(cache$pp, 2L, boundary, FUN = "<=")
  stopped <- rowSums(below) > 0L
  first_stop <- max.col(below, ties.method = "first")
  n_used <- rep(n_max, length(stopped))
  n_used[stopped] <- interim_counts[first_stop[stopped]]
  success <- !stopped & cache$final > pp_threshold
  success_n <- sum(success)
  stopped_n <- sum(stopped)
  c(
    success_n = success_n,
    rate = success_n / cache$attempted,
    pet = stopped_n / cache$attempted,
    ess = (sum(n_used) + cache$failed * n_max) / cache$attempted
  )
}

find_exact_target_pp <- function(cache, pp_grid, lambda_prime, interim_counts, n_max,
                                 target_success = 500L) {
  memo <- new.env(parent = emptyenv())
  count_at <- function(i) {
    key <- as.character(i)
    if (!exists(key, envir = memo, inherits = FALSE)) {
      value <- scenario_metrics(cache, pp_grid[i], lambda_prime, interim_counts, n_max)[["success_n"]]
      assign(key, as.integer(value), envir = memo)
    }
    get(key, envir = memo, inherits = FALSE)
  }
  n <- length(pp_grid)
  if (count_at(1L) < target_success || count_at(n) > target_success) return(integer(0L))

  lo <- 1L
  hi <- n
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (count_at(mid) <= target_success) hi <- mid else lo <- mid + 1L
  }
  first <- lo
  if (count_at(first) != target_success) return(integer(0L))

  lo <- first
  hi <- n
  while (lo < hi) {
    mid <- ceiling((lo + hi) / 2)
    if (count_at(mid) >= target_success) lo <- mid else hi <- mid - 1L
  }
  last <- lo
  if (count_at(last) != target_success) return(integer(0L))
  seq.int(first, last)
}

screen_one <- function(nm) {
  info <- project_info[project_info$project == nm, , drop = FALSE]
  if (nrow(info) != 1L) stop("Missing project metadata: ", nm)
  method <- info$method
  case <- info$case
  n_max <- info$n_max
  interim_counts <- case_interims[[case]]
  h0_path <- cache_table$path[cache_table$project == nm & cache_table$scenario == "H0"]
  h1_path <- cache_table$path[cache_table$project == nm & cache_table$scenario == "H1_standard"]
  if (length(h0_path) != 1L || length(h1_path) != 1L) stop("Cache mapping failed: ", nm)
  h0 <- prepare_cache(h0_path)
  h1 <- prepare_cache(h1_path)

  if (method == "RMUT") {
    pp_grid <- round(seq(0.80, 0.99, by = 0.0005), 4L)
    lambda_grid <- round(seq(0, 5, by = 0.005), 4L)
  } else if (method == "RMST") {
    pp_grid <- round(seq(0.80, 0.99, by = 0.001), 4L)
    lambda_grid <- round(seq(0, 15, by = 0.01), 4L)
  } else {
    pp_grid <- round(seq(0.80, 0.99, by = 0.001), 4L)
    lambda_grid <- round(seq(0, 5, by = 0.005), 4L)
  }

  rows <- vector("list", length(lambda_grid))
  for (j in seq_along(lambda_grid)) {
    idx <- find_exact_target_pp(h0, pp_grid, lambda_grid[j], interim_counts, n_max, 500L)
    if (!length(idx)) next
    rows[[j]] <- do.call(rbind, lapply(idx, function(i) {
      m0 <- scenario_metrics(h0, pp_grid[i], lambda_grid[j], interim_counts, n_max)
      m1 <- scenario_metrics(h1, pp_grid[i], lambda_grid[j], interim_counts, n_max)
      data.frame(
        project = nm,
        case = case,
        method = method,
        candidate_id = i + (j - 1L) * length(pp_grid),
        PP_final = pp_grid[i],
        lambda_prime = lambda_grid[j],
        H0_T1E = unname(m0["rate"]),
        H0_PET = unname(m0["pet"]),
        H0_ESS = unname(m0["ess"]),
        H1_Power = unname(m1["rate"]),
        H1_PET = unname(m1["pet"]),
        H1_ESS = unname(m1["ess"]),
        H0_go_n = as.integer(m0["success_n"]),
        H1_go_n = as.integer(m1["success_n"]),
        N_attempted_H0 = h0$attempted,
        N_attempted_H1 = h1$attempted,
        N_valid_H0 = h0$valid,
        N_valid_H1 = h1$valid,
        stringsAsFactors = FALSE
      )
    }))
  }
  exact <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(exact) || !nrow(exact)) stop("No exact-5% candidate on fine grid: ", nm)
  exact <- exact[exact$H1_Power >= 0.80, , drop = FALSE]
  if (!nrow(exact)) stop("No exact-5% candidate has H1 power >=80%: ", nm)
  exact[order(exact$PP_final, exact$lambda_prime, exact$candidate_id), , drop = FALSE]
}


prepare_cache_full <- function(path) {
  x <- readRDS(path)
  attempted <- length(x)
  if (attempted != 10000L) stop("Expected 10,000 records: ", path)
  valid <- vapply(x, function(z) {
    flag <- if (!is.null(z$ok)) isTRUE(z$ok) else if (!is.null(z$fit_ok)) isTRUE(z$fit_ok) else TRUE
    pp <- if (!is.null(z$PP)) z$PP else z$PP_interims
    flag && length(pp) == 3L && all(is.finite(pp)) &&
      length(z$PP_final) == 1L && is.finite(z$PP_final)
  }, logical(1L))
  if (mean(valid) < 0.95) stop("Valid-fit rate below 95%: ", path)
  pp <- matrix(NA_real_, nrow = attempted, ncol = 3L)
  valid_index <- which(valid)
  pp[valid_index, ] <- do.call(rbind, lapply(x[valid], function(z) {
    as.numeric(if (!is.null(z$PP)) z$PP else z$PP_interims)
  }))
  list(pp = pp, valid = valid, attempted = attempted, path = path)
}

n_used_vector <- function(cache, pp_threshold, lambda_prime, interim_counts, n_max) {
  n_used <- rep(as.numeric(n_max), cache$attempted)
  valid_index <- which(cache$valid)
  boundary <- pp_threshold * (interim_counts / n_max)^lambda_prime
  below <- sweep(cache$pp[valid_index, , drop = FALSE], 2L, boundary, FUN = "<=")
  stopped <- rowSums(below) > 0L
  if (any(stopped)) {
    first_stop <- max.col(below[stopped, , drop = FALSE], ties.method = "first")
    n_used[valid_index[stopped]] <- interim_counts[first_stop]
  }
  n_used
}

