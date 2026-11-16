# Portable input/output helpers for Cox/RMST parameter selection.

selector_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg) != 1L) stop("Run the entry with Rscript.")
  dirname(normalizePath(sub("^--file=", "", arg), mustWork = TRUE, winslash = "/"))
}

selector_args <- function(usage, require_eligible = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || "--help" %in% args) {
    cat(usage, "\n")
    quit(save = "no", status = if (length(args)) 0L else 1L)
  }
  result <- list()
  i <- 1L
  allowed <- c("cache-manifest", "output-dir", "eligible", "case", "method", "reference")
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key == "verify-reference") {
      if (!is.null(result[[key]])) stop("Duplicate argument: ", args[i])
      result[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(args[i], "--") || !key %in% allowed ||
        i == length(args) || startsWith(args[i + 1L], "--")) {
      stop("Unknown or incomplete argument: ", args[i], "\n", usage)
    }
    if (!is.null(result[[key]])) stop("Duplicate argument: ", args[i])
    result[[key]] <- args[i + 1L]
    i <- i + 2L
  }
  required <- c("cache-manifest", "output-dir", if (require_eligible) "eligible")
  if (any(!required %in% names(result))) stop("Required arguments: ", paste(required, collapse = ", "))
  result
}

selector_case <- function(x) {
  x <- as.character(x)
  x[x == "IMpower"] <- "IMpower131"
  x[x == "KEYNOTE598"] <- "KEYNOTE-598"
  x[x == "OAK"] <- "OAK-v5"
  x
}

selector_manifest <- function(path, method = NULL, case = NULL) {
  path <- normalizePath(path, mustWork = TRUE, winslash = "/")
  tab <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("method", "case", "scenario", "path")
  if (!all(required %in% names(tab)) || !nrow(tab)) stop("Manifest requires nonempty method,case,scenario,path columns.")
  if (anyNA(tab[, required]) || any(!nzchar(tab$path))) stop("Missing manifest values.")
  tab$case <- selector_case(tab$case)
  tab$scenario[tab$scenario %in% c("H1-standard", "H1_calibration")] <- "H1_standard"
  cases <- c("IMpower131", "KEYNOTE-598", "OAK-v5", "survival_colon")
  if (any(!tab$method %in% c("Cox", "RMST")) || any(!tab$case %in% cases) ||
      any(!tab$scenario %in% c("H0", "H1_standard"))) {
    stop("Manifest accepts Cox/RMST, the four formal cases, and H0/H1_standard only; H1-observed is excluded.")
  }
  if (!is.null(method)) {
    if (!method %in% c("Cox", "RMST")) stop("--method must be Cox or RMST.")
    tab <- tab[tab$method == method, , drop = FALSE]
  }
  if (!is.null(case)) {
    case <- selector_case(case)
    if (!case %in% cases) stop("Unknown --case.")
    tab <- tab[tab$case == case, , drop = FALSE]
  }
  if (!nrow(tab)) stop("No selected manifest rows.")
  key <- paste(tab$method, tab$case, tab$scenario, sep = "/")
  if (anyDuplicated(key)) stop("Duplicate method/case/scenario rows in manifest.")
  groups <- unique(tab[, c("method", "case"), drop = FALSE])
  for (g in seq_len(nrow(groups))) {
    z <- tab[tab$method == groups$method[g] & tab$case == groups$case[g], , drop = FALSE]
    if (!setequal(z$scenario, c("H0", "H1_standard")) || nrow(z) != 2L) {
      stop("Each selected method/case must have one H0 and one H1_standard cache.")
    }
  }
  resolve <- function(p) {
    p <- path.expand(p)
    if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", p)) p <- file.path(dirname(path), p)
    normalizePath(p, mustWork = TRUE, winslash = "/")
  }
  tab$path <- vapply(tab$path, resolve, character(1L))
  tab$project <- paste(tab$method, gsub("-", "", tab$case, fixed = TRUE), sep = "_")
  # Names are identifiers only; the historical case labels retain OAK-v5.
  groups$project <- paste(groups$method, gsub("-", "", groups$case, fixed = TRUE), sep = "_")
  groups$n_max <- c(IMpower131 = 680L, `KEYNOTE-598` = 568L, `OAK-v5` = 850L, survival_colon = 620L)[groups$case]
  interims <- list(IMpower131 = c(272L,408L,544L), `KEYNOTE-598` = c(228L,340L,454L),
                   `OAK-v5` = c(340L,510L,680L), survival_colon = c(248L,372L,496L))
  list(project_info = groups, cache_table = tab[, c("project","scenario","path")],
       case_interims = interims, manifest_path = path)
}

selector_new_output <- function(path, package_root) {
  source(file.path(package_root,"src","common","output_paths.R"),local=TRUE)
  path <- external_output_path(path,package_root)
  package_root <- normalizePath(package_root, mustWork = TRUE, winslash = "/")
  if (identical(tolower(path), tolower(package_root)) ||
      startsWith(tolower(path), paste0(tolower(package_root), "/"))) {
    stop("Output directory must be outside the source package.")
  }
  if (file.exists(path) || dir.exists(path)) stop("Use a new output directory: ", path)
  if (!dir.create(path, recursive = FALSE)) stop("Cannot create output directory: ", path)
  path
}

selector_check_reference <- function(selected, reference) {
  reference <- utils::read.csv(normalizePath(reference, mustWork = TRUE), stringsAsFactors = FALSE)
  if (!all(c("method","case","candidate_id","c","gamma") %in% names(reference))) stop("Invalid reference table.")
  reference$case <- selector_case(reference$case)
  for (i in seq_len(nrow(selected))) {
    z <- selected[i, , drop = FALSE]
    r <- reference[reference$method == z$method & reference$case == z$case, , drop = FALSE]
    if (nrow(r) != 1L || r$candidate_id != z$candidate_id ||
        abs(r$c - z$PP_final) > 1e-12 || abs(r$gamma - z$lambda_prime) > 1e-12) {
      stop("Computed selection differs from the locked reference for ", z$case, " / ", z$method,
           ". No reference parameters were substituted.")
    }
  }
  invisible(TRUE)
}
