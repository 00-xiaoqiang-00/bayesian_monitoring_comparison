# Resolve an output path without creating it, then enforce the package boundary.
external_output_path <- function(path, code_root) {
  if (!nzchar(path)) stop("Output destination must not be empty.")
  requested <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\]|\\\\)", requested)) requested <- file.path(getwd(), requested)
  ancestor <- requested
  suffix <- character()
  while (!file.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) stop("Cannot resolve output destination.")
    suffix <- c(basename(ancestor), suffix)
    ancestor <- parent
  }
  if (!dir.exists(ancestor)) stop("An output ancestor is not a directory.")
  resolved <- normalizePath(ancestor, winslash="/", mustWork=TRUE)
  for (part in suffix) {
    if (part == ".") next
    resolved <- if (part == "..") dirname(resolved) else file.path(resolved, part)
    if (file.exists(resolved)) resolved <- normalizePath(resolved,winslash="/",mustWork=TRUE)
  }
  code_root <- normalizePath(code_root,winslash="/",mustWork=TRUE)
  root <- if (dir.exists(file.path(dirname(code_root),"04_Data_and_Results"))) dirname(code_root) else code_root
  key <- sub("/+$", "", gsub("\\", "/", resolved, fixed=TRUE))
  boundary <- sub("/+$", "", gsub("\\", "/", root, fixed=TRUE))
  if (.Platform$OS.type == "windows") { key <- tolower(key); boundary <- tolower(boundary) }
  if (identical(key,boundary) || startsWith(key,paste0(boundary,"/"))) stop("Choose an output destination outside the submission/code package.")
  resolved
}
