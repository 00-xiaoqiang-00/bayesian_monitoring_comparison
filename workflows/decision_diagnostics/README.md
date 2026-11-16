# OAK posterior-probability diagnostic (Table 6)

`summarize_oak_table6.R` reconstructs Table 6 for OAK under H1-observed from the saved four-look posterior-probability (PP) caches of the four locked independent-validation designs. It uses base R and performs no model fitting or patient/trial simulation.

## Inputs and execution

Run from the repository root, with the four generated RDS caches in an external directory:

```text
Rscript workflows/decision_diagnostics/summarize_oak_table6.R --cache-dir /path/to/saved_pp_caches --output-dir /path/to/new_external_results
```

The default input tables are `config/current_parameters.csv`, `config/trial_process.csv`, and `config/analysis_schedule.csv`. The default cache filenames and formats are specified in [`oak_table6_cache_spec.csv`](../../config/oak_table6_cache_spec.csv). These are generated validation caches, not raw trial data. Their large RDS files are not bundled in this package. To use caches saved with other filenames, supply a CSV with the same three columns using `--cache-spec /path/to/cache_spec.csv`; its `relative_cache` paths are resolved against `--cache-dir`.

The cache formats are:

- `stats_list` (Cox, RMST, RMUT): a list of replicate records, each with three `PP` values and one `PP_final`; the script respects `ok` or `fit_ok` when present.
- `winratio_matrix`: an RDS list containing a four-column `PP` matrix, with optional `look_cohort` and `N.max` metadata checked against the locked trial schedule.

Generate the four caches with the formal validation entries listed in the root analysis index.

Explicit input overrides are available with `--parameters`, `--trial-process`, and `--schedule`. Run `--help` for the complete syntax. The output directory must be outside the repository. Existing output files are not overwritten.

## Definitions and outputs

The valid-only denominator requires a true fit flag when one is supplied and finite PP at all four looks. All four archived OAK caches have 10,000 valid replicates. The main interim information fractions are derived from the locked cohort schedule and checked against `trial_process.csv`.

At interim look k the futility rule is `PP_k <= eta * t_k^gamma`; subsequent stops are evaluated only among replicates continuing past earlier looks. Final go requires both reaching the final analysis and `PP_final > eta`. Equality at the final threshold is not success. The diagnostic final-success proportion uses `PP_final > eta` for every valid replicate, including those already stopped at an interim analysis. This diagnostic proportion is neither the go probability of the full locked monitoring strategy nor an operating characteristic of a design recalibrated without interim stopping.

Outputs are:

- `OAK_H1observed_table6.csv`: four rows with PP medians, futility boundaries, sequential stop counts, final counts, actual go/PET/ESS, diagnostic final-success counts/proportions, and threshold-tie counts. Probabilities are stored as fractions, not percentages.
- `OAK_H1observed_PP_quantiles.csv`: 5th/25th/50th/75th/95th PP percentiles at all four looks, using R's default type-7 quantiles.
- `input_provenance.csv`: resolved input paths and MD5 checksums, written only into the external output directory.

To compare a rerun against the archived four-row summary, add `--reference /path/to/OAK_H1observed_table6.csv`. The script compares every column and uses a numerical tolerance of 1e-12. The archived results are under `04_Data_and_Results/Decision_Diagnostics` in the full submission distribution.

The archived actual go counts are 2,687 (Win Ratio), 0 (Cox), 0 (RMST), and 6 (RMUT). The diagnostic final-success proportions are 0.9394, 0.1993, 0.5676, and 0.6873, respectively. No replicate reaching the final analysis equals its final threshold; therefore strict final-threshold handling does not change these actual go counts.
