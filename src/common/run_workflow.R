v33_code_root <- function() {
  configured <- Sys.getenv("V33_CODE_ROOT", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Run a workflow using Rscript.")
  directory <- dirname(normalizePath(sub("^--file=", "", hit), winslash = "/", mustWork = TRUE))
  repeat {
    if (file.exists(file.path(directory, "config", "current_parameters.csv"))) return(directory)
    parent <- dirname(directory)
    if (identical(parent, directory)) stop("Cannot locate the code-package root.")
    directory <- parent
  }
}

v33_run <- function(implementation, method, case, analysis, fixed_arguments = character()) {
  if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_ALL", "English_United States.utf8")))
  root <- v33_code_root()
  args <- commandArgs(TRUE)
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript WORKFLOW.R --output-dir EXTERNAL_DIRECTORY [analysis options]\n",
        "       Rscript WORKFLOW.R --preflight-only\n",
        "Win Ratio calibration/validation/interim workflows accept --prepare-only with an external output directory to copy dependency sources without compilation or analysis.\n",
        "The preflight checks supplied source paths, R syntax and fixed parameters; it performs no model fitting.\n",
        "Dependence runs for Cox/RMST/RMUT require --level=low|high.\n",
        "Win Ratio dependence runs also require --runtime-dir=PREPARED_RUNTIME, --case=CASE, --level=LEVEL and --mode=formal.\n", sep = "")
    quit(status = 0L)
  }
  option <- function(name, default = "") {
    prefix <- paste0("--", name, "=")
    hit <- args[startsWith(args, prefix)]
    separate <- which(args == paste0("--", name))
    if (length(hit) + length(separate) > 1L) stop("Duplicate option: ", name)
    if (length(hit)) return(substring(hit, nchar(prefix) + 1L))
    if (length(separate)) {
      if (separate == length(args)) stop("Missing value for ", name)
      return(args[separate + 1L])
    }
    default
  }
  remove_option <- function(name) {
    separate <- which(args == paste0("--", name))
    remove <- c(which(startsWith(args, paste0("--", name, "="))), separate,
                if (length(separate)) separate + 1L else integer())
    if (length(remove)) args <<- args[-remove]
  }
  output <- option("output-dir", Sys.getenv("V33_OUTPUT_DIR", unset = ""))
  runtime <- option("runtime-dir")
  prepare_only <- "--prepare-only" %in% args
  args <- args[args != "--prepare-only"]
  level <- option("level", Sys.getenv("SENSITIVITY_RHO_LEVEL", unset = ""))
  remove_option("output-dir"); remove_option("runtime-dir")
  if (method != "Win Ratio") remove_option("level")
  target <- file.path(root, implementation)
  if (!file.exists(target)) stop("Missing method implementation: ", target)
  required <- c("config/current_parameters.csv", "config/dgm_parameters.csv", "config/trial_process.csv",
                "config/analysis_schedule.csv", "config/rmut_parameters.tsv", "config/rho_parameters.csv")
  if (any(!file.exists(file.path(root, required)))) stop("Missing fixed input: ", paste(required[!file.exists(file.path(root, required))], collapse = ", "))
  parameters <- read.csv(file.path(root, "config/current_parameters.csv"), stringsAsFactors = FALSE)
  if (any(!is.finite(parameters$c)) || any(parameters$c <= 0 | parameters$c >= 1) || any(!is.finite(parameters$gamma)) || any(parameters$gamma < 0)) stop("Invalid locked parameters.")
  normalize_case <- function(value) {
    value[value == "IMpower"] <- "IMpower131"
    value[value == "KEYNOTE598"] <- "KEYNOTE-598"
    value
  }
  if (method != "Win Ratio" && analysis != "calibration") {
    row <- parameters[parameters$method == method & normalize_case(parameters$case) == normalize_case(case), , drop = FALSE]
    if (nrow(row) != 1L) stop("Expected one locked method/case parameter row.")
    text <- paste(readLines(target, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    scalar <- function(name) {
      hit <- regmatches(text, regexec(paste0(name, "[[:space:]]*<-[[:space:]]*([.0-9]+)"), text))[[1]]
      if (length(hit) != 2L) stop("Cannot verify the locked parameter in the implementation: ", name)
      as.numeric(hit[2L])
    }
    if (abs(scalar("iv_locked_pp") - row$c) > 1e-12 || abs(scalar("iv_locked_lambda") - row$gamma) > 1e-12) stop("Implementation and locked configuration thresholds disagree.")
  }
  if (method == "Win Ratio") {
    wr <- read.delim(file.path(root, "config/win_ratio_parameters.tsv"), stringsAsFactors = FALSE)
    reference <- parameters[parameters$method == "Win Ratio", , drop = FALSE]
    matched <- match(normalize_case(wr$case_key), normalize_case(reference$case))
    if (nrow(wr) != 4L || anyNA(matched) || any(abs(wr$c_F - reference$c[matched]) > 1e-12) || any(abs(wr$gamma_E - reference$gamma[matched]) > 1e-12)) stop("Win Ratio locked configurations disagree.")
  }
  parse(target, encoding = "UTF-8")
  if ("--preflight-only" %in% args) {
    scope <- if (method == "Win Ratio") "src/methods/win_ratio" else file.path("src/methods", tolower(method))
    files <- unique(c(target, list.files(file.path(root, scope), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
                      list.files(file.path(root, "src/common"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE),
                      list.files(file.path(root, "src/dgm"), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)))
    for (file in files) parse(file, encoding = "UTF-8")
    if (method == "Win Ratio") {
      if (!file.exists(file.path(root, "config/win_ratio_parameters.tsv"))) stop("Missing Win Ratio locked parameters.")
      vendor <- file.path(root, "third_party/win_ratio")
      if (!dir.exists(vendor)) stop("Missing third-party source.")
    }
    cat("INPUT_PREFLIGHT_OK method=", method, " case=", case, " analysis=", analysis,
        " parsed_sources=", length(files), "; no model loading, simulation or MCMC\n", sep = "")
    quit(status = 0L)
  }
  if (!nzchar(output)) stop("Supply --output-dir outside the submission package.")
  output <- path.expand(output)
  parent <- normalizePath(dirname(output), winslash = "/", mustWork = TRUE)
  output <- file.path(parent, basename(output))
  protected <- if (dir.exists(file.path(dirname(root), "04_Data_and_Results"))) dirname(root) else root
  canonical <- function(x) tolower(gsub("\\\\", "/", x))
  if (identical(canonical(output), canonical(protected)) || startsWith(paste0(canonical(output), "/"), paste0(canonical(protected), "/"))) stop("Output must be outside the entire submission.")
  Sys.setenv(V33_CODE_ROOT = root, V33_OUTPUT_DIR = output, PROJECT_ROOT = root,
             RMUT_LOCKED_PARAMETER_MANIFEST = file.path(root, "config/rmut_parameters.tsv"),
             RHO_PARAMETER_FILE = file.path(root, "config/rho_parameters.csv"),
             RHO_SENSITIVITY_HELPER = file.path(root, "src/common/rho_parameters.R"),
             RHO_SENSITIVITY_OUTPUT_ROOT = output, COX_COLON_CALIBRATION_OUTPUT_ROOT = output,
             WINRATIO_OUTPUT_ROOT = file.path(output, "case_setup"),
             WINRATIO_VALIDATION_OUTPUT_ROOT = output, WINRATIO_S8_OUTPUT_ROOT = output)
  case_key <- switch(case, IMpower = "IMpower131", `KEYNOTE-598` = "KEYNOTE598", case)
  source_case <- switch(case_key, IMpower131 = "IMpower", KEYNOTE598 = "KEYNOTE598", case_key)
  if (method == "Win Ratio" && analysis %in% c("calibration", "independent validation", "50/80 interim sensitivity")) {
    source_packages <- file.path(root, "third_party/win_ratio/packages")
    runtime_packages <- file.path(output, "runtime_packages")
    dir.create(runtime_packages, recursive = TRUE, showWarnings = FALSE)
    packages <- list.dirs(source_packages, recursive = FALSE, full.names = TRUE)
    for (package in packages) {
      destination <- file.path(runtime_packages, basename(package))
      if (!dir.exists(destination)) {
        if (!file.copy(package, runtime_packages, recursive = TRUE)) stop("Cannot prepare external package: ", package)
      } else {
        files <- list.files(package, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
        files <- files[!file.info(files)$isdir]
        relative <- substring(files, nchar(package) + 2L)
        # pkgload regenerates these two interfaces in the external copy.
        keep <- !gsub("\\\\", "/", relative) %in% c("R/RcppExports.R", "src/RcppExports.cpp")
        targets <- file.path(destination, relative[keep])
        if (any(!file.exists(targets)) || !identical(unname(tools::md5sum(files[keep])), unname(tools::md5sum(targets)))) stop("External package source differs; choose a new output directory.")
      }
    }
    Sys.setenv(WINRATIO_PACKAGE_ROOT = runtime_packages)
    if (case != "all four") Sys.setenv(BWR_REPO_ROOT = file.path(runtime_packages, paste0("winratio_", source_case)))
    if (prepare_only) {
      cat("DEPENDENCY_PREPARATION_OK external_packages=", length(packages), "; no compilation, model loading, simulation or MCMC\n", sep = "")
      quit(status = 0L)
    }
  } else if (prepare_only) stop("--prepare-only is available for Win Ratio calibration, validation and interim-time workflows.")
  if (analysis == "dependence sensitivity" && method != "Win Ratio") {
    if (!level %in% c("low", "high")) stop("Use --level=low or --level=high.")
    rho <- switch(case_key, IMpower131 = c(low = .3, high = .7), KEYNOTE598 = c(low = .3, high = .7), OAK = c(low = .5, high = .9), survival_colon = c(low = .7, high = .98))[[level]]
    if (is.null(rho)) stop("Unknown dependence case.")
    Sys.setenv(SENSITIVITY_RHO_LEVEL = level, SENSITIVITY_RHO = as.character(rho), RHO_SENSITIVITY_VERSION_DIR = file.path(root, "config"))
    if (method %in% c("Cox", "RMST")) Sys.setenv(RHO_SENSITIVITY_DATA_GENERATOR = file.path(root, "src/dgm/dependence", paste0(if (method == "Cox") "cox_" else "RMST_", source_case, ".R")))
  }
  if (method == "Win Ratio" && analysis == "dependence sensitivity") {
    if (!nzchar(runtime)) stop("Prepare an external Win Ratio runtime, then pass --runtime-dir=PATH.")
    target <- file.path(normalizePath(runtime, winslash = "/", mustWork = TRUE), "validate_winratio_rho.R")
    if (!file.exists(target)) stop("Prepared runtime lacks validate_winratio_rho.R.")
  }
  if (method == "Win Ratio" && analysis == "50/80 interim sensitivity") args <- c(args, paste0("--output-root=", output))
  if (analysis == "Win Ratio runtime preparation") fixed_arguments <- c(fixed_arguments, paste0("--output=", output))
  args <- c(fixed_arguments, args)
  executable <- file.path(R.home("bin"), "Rscript")
  status <- system2(executable, c("--vanilla", shQuote(target), vapply(args, shQuote, character(1))))
  quit(status = as.integer(status))
}
