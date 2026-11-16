#!/usr/bin/env Rscript
# Formal paired-MCSE final selection.
options(stringsAsFactors = FALSE)
entry_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(entry_arg) != 1L) stop("Run with Rscript.")
entry_dir <- dirname(normalizePath(sub("^--file=", "", entry_arg), mustWork = TRUE))
source(file.path(entry_dir, "..", "..", "..", "src", "common", "parameter_selection", "selection_io.R"), local = TRUE, encoding = "UTF-8")

main <- function() {
  a <- selector_args("Rscript select_cox_rmst.R --cache-manifest CACHE_MAP.tsv --eligible ELIGIBLE.tsv --output-dir NEW_EXTERNAL_DIRECTORY [--method Cox|RMST] [--case CASE] [--verify-reference] [--reference CSV]", require_eligible = TRUE)
  mapping <- selector_manifest(a[["cache-manifest"]], a[["method"]], a[["case"]])
  project_map <- mapping$project_info
  cache_table <- mapping$cache_table
  case_interims <- mapping$case_interims
  eligible <- utils::read.delim(normalizePath(a[["eligible"]], mustWork = TRUE), check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("case","method","candidate_id","PP_final","lambda_prime","H0_T1E","H1_Power","H0_PET","H0_ESS","H1_PET","H1_ESS")
  if (!all(required %in% names(eligible))) stop("Eligible table lacks required metric columns.")
  if (any(!eligible$method %in% c("Cox","RMST"))) eligible <- eligible[eligible$method %in% c("Cox","RMST"),,drop=FALSE]
  eligible$case <- selector_case(eligible$case)
  group_key <- paste(project_map$method, project_map$case, sep = "/")
  eligible <- eligible[paste(eligible$method, eligible$case, sep = "/") %in% group_key &
                         eligible$H0_T1E == 0.05 & eligible$H1_Power >= 0.80 &
                         eligible$PP_final >= 0.85 & eligible$lambda_prime >= 0.01, , drop = FALSE]
  if (!nrow(eligible)) stop("No candidates pass the exact-T1E, power, c and gamma gates.")
  if (anyDuplicated(paste(eligible$method, eligible$case, eligible$candidate_id, sep = "/"))) stop("Duplicate candidate IDs within a group.")
  source(file.path(entry_dir, "..", "..", "..", "src", "common", "parameter_selection", "selection_functions.R"), local = TRUE, encoding = "UTF-8")
  output_dir <- selector_new_output(a[["output-dir"]], file.path(entry_dir, "..", "..", ".."))
  case_order <- c("IMpower131", "KEYNOTE-598", "OAK-v5", "survival_colon")
  method_order <- c("Cox", "RMST")
all_results <- vector("list", nrow(project_map))
summary_rows <- vector("list", nrow(project_map))
selected_rows <- vector("list", nrow(project_map))

for (g in seq_len(nrow(project_map))) {
  info <- project_map[g, , drop = FALSE]
  key_rows <- eligible$case == info$case & eligible$method == info$method
  candidates <- eligible[key_rows, , drop = FALSE]
  if (!nrow(candidates)) stop("No candidates for ", info$case, " / ", info$method)

  h0_path <- cache_table$path[
    cache_table$project == info$project & cache_table$scenario == "H0"
  ]
  if (length(h0_path) != 1L) stop("H0 cache mapping failed for ", info$project)
  cache <- prepare_cache_full(h0_path)
  interim_counts <- case_interims[[info$case]]
  n_max <- info$n_max

  strict_order <- order(
    candidates$H0_ESS,
    -candidates$PP_final,
    -candidates$H1_Power,
    -candidates$H0_PET,
    candidates$candidate_id
  )
  candidates <- candidates[strict_order, , drop = FALSE]
  strict_min <- candidates[1L, , drop = FALSE]
  n_min <- n_used_vector(
    cache, strict_min$PP_final, strict_min$lambda_prime,
    interim_counts, n_max
  )
  min_recomputed <- mean(n_min)
  if (abs(min_recomputed - strict_min$H0_ESS) > 1e-8) {
    stop("Strict-min H0 ESS mismatch for ", info$case, " / ", info$method,
         ": stored=", strict_min$H0_ESS, ", recomputed=", min_recomputed)
  }

  recomputed <- numeric(nrow(candidates))
  delta <- numeric(nrow(candidates))
  paired_se <- numeric(nrow(candidates))
  near <- logical(nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    n_i <- n_used_vector(
      cache, candidates$PP_final[i], candidates$lambda_prime[i],
      interim_counts, n_max
    )
    recomputed[i] <- mean(n_i)
    if (abs(recomputed[i] - candidates$H0_ESS[i]) > 1e-8) {
      stop("H0 ESS mismatch for candidate ", candidates$candidate_id[i],
           " in ", info$case, " / ", info$method,
           ": stored=", candidates$H0_ESS[i], ", recomputed=", recomputed[i])
    }
    diff <- n_i - n_min
    delta[i] <- mean(diff)
    paired_se[i] <- stats::sd(diff) / sqrt(length(diff))
    near[i] <- delta[i] <= 2 * paired_se[i] + 1e-12
  }

  candidates$strict_min_candidate_id <- strict_min$candidate_id
  candidates$H0_ESS_recomputed <- recomputed
  candidates$delta_H0_ESS_vs_min <- delta
  candidates$paired_SE_delta <- paired_se
  candidates$near_optimal_2SE <- near
  candidates$boundary_flag <- (
    candidates$PP_final <= 0.855 |
      candidates$lambda_prime <= 0.015
  )

  near_set <- candidates[candidates$near_optimal_2SE, , drop = FALSE]
  if (!nrow(near_set)) stop("Near-optimal set is empty")
  robust_order <- order(
    -near_set$H1_Power,
    -near_set$PP_final,
    near_set$H1_PET,
    near_set$boundary_flag,
    near_set$candidate_id
  )
  near_set <- near_set[robust_order, , drop = FALSE]
  chosen <- near_set[1L, , drop = FALSE]
  candidates$selected_robust <- candidates$candidate_id == chosen$candidate_id
  candidates$selection_rank_within_near <- NA_integer_
  candidates$selection_rank_within_near[
    match(near_set$candidate_id, candidates$candidate_id)
  ] <- seq_len(nrow(near_set))

  changed <- chosen$candidate_id != strict_min$candidate_id
  summary_rows[[g]] <- data.frame(
    case = info$case,
    method = info$method,
    n_eligible = nrow(candidates),
    n_near_optimal_2SE = nrow(near_set),
    strict_min_candidate_id = strict_min$candidate_id,
    strict_min_PP_final = strict_min$PP_final,
    strict_min_lambda_prime = strict_min$lambda_prime,
    strict_min_H1_Power = strict_min$H1_Power,
    strict_min_H0_ESS = strict_min$H0_ESS,
    selected_candidate_id = chosen$candidate_id,
    selected_PP_final = chosen$PP_final,
    selected_lambda_prime = chosen$lambda_prime,
    selected_H1_Power = chosen$H1_Power,
    selected_H0_PET = chosen$H0_PET,
    selected_H0_ESS = chosen$H0_ESS,
    selected_H1_PET = chosen$H1_PET,
    selected_H1_ESS = chosen$H1_ESS,
    selected_delta_H0_ESS = chosen$delta_H0_ESS_vs_min,
    selected_paired_SE_delta = chosen$paired_SE_delta,
    changed_from_strict_min = changed,
    selected_boundary_flag = chosen$boundary_flag,
    stringsAsFactors = FALSE
  )
  chosen$selection_status <- "near_optimal_2SE_then_robust_ranking"
  selected_rows[[g]] <- chosen
  candidates$project <- info$project
  candidates$H0_cache_path <- h0_path
  all_results[[g]] <- candidates
  message(
    info$case, " / ", info$method,
    ": eligible=", nrow(candidates),
    ", near=", nrow(near_set),
    ", strict=", strict_min$candidate_id,
    ", selected=", chosen$candidate_id,
    ", changed=", changed
  )
}

all_results <- do.call(rbind, all_results)
summary_table <- do.call(rbind, summary_rows)
selected <- do.call(rbind, selected_rows)


  if (isTRUE(a[["verify-reference"]])) {
    reference <- a[["reference"]]
    if (is.null(reference)) reference <- file.path(entry_dir, "..", "..", "..", "config", "current_parameters.csv")
    selector_check_reference(selected, reference)
  }
  summary_table <- summary_table[order(match(summary_table$case, case_order), match(summary_table$method, method_order)),,drop=FALSE]
  selected <- selected[order(match(selected$case, case_order), match(selected$method, method_order)),,drop=FALSE]
  all_results <- all_results[order(match(all_results$case, case_order), match(all_results$method, method_order), all_results$selection_rank_within_near, all_results$H0_ESS),,drop=FALSE]
  utils::write.table(summary_table, file.path(output_dir, "near_optimal_group_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(selected, file.path(output_dir, "selected_parameter_groups.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(all_results, file.path(output_dir, "all_candidates_with_paired_mcse.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(c(
    "algorithm=NEAR_OPTIMAL_PAIRED_MCSE_V1_20260806;scope=Cox_RMST",
    "H0_T1E=0.05_exact;H1_power_min=0.80;c_min=0.85;gamma_min=0.01",
    "near_optimal_rule=delta_H0_ESS_le_2_times_paired_SE_delta",
    "ranking=H1_Power_desc,c_desc,H1_PET_asc,boundary_flag_asc,candidate_id_asc",
    "boundary_flag=c_le_0.855_or_gamma_le_0.015",
    "H1_observed_used=false;locked_parameters_used_for_selection=false",
    paste0("reference_compared_after_selection=", isTRUE(a[["verify-reference"]])),
    paste0("input_manifest=", mapping$manifest_path),
    paste0("input_manifest_md5=", unname(tools::md5sum(mapping$manifest_path)))), file.path(output_dir, "SELECTION_METADATA.txt"))
  cat("SELECTION_COMPLETE projects=", nrow(selected), " changed_from_strict_min=", sum(summary_table$changed_from_strict_min), "\n", sep = "")
}
if (sys.nframe() == 0L) main()
