# Cox/RMST formal parameter selection

These entries implement exact-count calibration screening and paired-MCSE parameter selection for Cox and RMST. Inputs are supplied through command-line arguments, and generated files are written to external output directories.

The selection is a two-stage process:

1. `screen_cox_rmst.R` reads H0 and H1-standard posterior-probability caches from a tab-delimited manifest. It checks 10,000 attempted records, requires a valid-fit rate of at least 95%, finds empirical H0 success count exactly 500/10,000 on the method-specific fine grid, and keeps H1-standard power at least 80%.
2. `select_cox_rmst.R` filters the eligible table using H0 T1E exactly 0.05, H1-standard power at least 80%, `c` at least 0.85 and `gamma` at least 0.01. It recomputes H0 ESS from the same H0 cache, defines the paired-MCSE near-optimal set by `delta_H0_ESS <= 2 * paired_SE_delta`, and ranks that set by H1-standard power descending, final cutoff `c` (`PP_final`) descending, H1-standard PET ascending, diagnostic boundary flag ascending, and candidate ID ascending. The boundary flag is true when `c <= 0.855` or `gamma <= 0.015`. The optional `--verify-reference` check compares the selected rows with `config/current_parameters.csv`; it does not substitute locked values for a computed result.

The manifest must contain `method`, `case`, `scenario`, and `path` columns, with one `H0` and one `H1_standard` row for each selected method/case. Relative cache paths resolve relative to the manifest file. `H1_observed` is excluded from selection.

Example for one case, using caches outside the package:

```text
Rscript workflows/parameter_selection/cox_rmst/screen_cox_rmst.R \
  --cache-manifest /path/to/cache_manifest.tsv \
  --output-dir /path/to/new/screen_output \
  --method Cox --case survival_colon

Rscript workflows/parameter_selection/cox_rmst/select_cox_rmst.R \
  --cache-manifest /path/to/cache_manifest.tsv \
  --eligible /path/to/new/screen_output/eligible_parameter_groups.tsv \
  --output-dir /path/to/new/selection_output \
  --method Cox --case survival_colon --verify-reference
```

Use new output directories outside the source package. The selector reads saved PP caches; it does not fit models, generate patients, run MCMC, or run independent validation. The cache-generation entries are under `workflows/calibration/` and are listed in the root `analysis_index.csv`.

`cache_manifest_template.tsv` is a layout template only and contains no result cache. The computed selections can be checked against the current eight Cox/RMST reference rows using `--verify-reference`.
