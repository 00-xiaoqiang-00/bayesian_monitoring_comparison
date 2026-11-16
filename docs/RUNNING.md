# Running analyses

Run commands from the repository root. Use
`analysis_index.csv` to choose the workflow for the desired method, case and
analysis. Configurations are under `config/`; fixed inputs are under `inputs/`.
Archived result tables are not bundled here. If a verification mode requires them, supply the results archive outside this repository and pass its documented path.

## Dependencies

Python workflows require Python 3.10 or later and the packages in the root
`requirements.txt`:

```text
python -m pip install -r requirements.txt
```

The calculation configurations record numerical methods and software-version
information where relevant. Monte Carlo streams and numerical libraries can
affect bitwise reproduction; use the recorded settings for numerical work.

R workflows use R 4.4.x. Install packages required by the selected workflow:
`brms`, `rstan`, `posterior`, `survival`, `MASS`, `Rcpp`, `RcppArmadillo`,
`RcppDist`, `mvtnorm`, `numDeriv`, `pkgload`, `digest` and `jsonlite`, as applicable. Stan and included C++
sources require a working compiler toolchain. The individual workflow README
identifies any additional fitting or reconstruction requirements.

## Execution modes

- Population-effect archive verification checks the saved result files,
  configurations and numerical definitions without new sampling.
- Explicit population-effect recomputation performs the documented numerical
  integration or batch Monte Carlo calculation. The optional integrals-only
  mode computes RMST/RMUT without WR/Cox sampling.
- Frozen-DGM assessment evaluates supplied parameters against fixed curves.
  DGM fitting and marginal reconstruction are separate workflows.
- Parameter selection and PP diagnostics read saved PP paths and do not fit
  models or simulate patients.
- Calibration, independent validation and sensitivity workflows run the
  corresponding monitoring analyses with their locked settings.

## Monitoring workflows

Public R monitoring entries accept `--preflight-only` for an input/source
check and require an external output directory for an actual run:

```text
Rscript workflows/validation/run_cox_oak.R --preflight-only
Rscript workflows/validation/run_cox_oak.R --output-dir /path/to/new_external_results
```

Dependence workflows additionally select `--level=low` or `--level=high`.
Win Ratio dependence runs use an externally prepared runtime. Read the
monitoring workflow instructions before launching a formal run; the original
replicate counts, seed streams, fitting settings and retry rules are part of
the analysis definition.

Keep all generated output outside the entire repository. PP caches,
compiled runtimes and model fits needed by a subsequent analysis are generated
there and supplied through that workflow's documented input arguments.

RMUT ESS post-processing expects the supplied data directory to contain
`RMUT_Run_Summaries.csv` directly. Write the derived ESS view outside the
repository:

```text
python -B workflows/postprocessing/summarize_rmut_ess_views.py /path/to/external_summary_directory /path/to/new_external_ess_views.csv
```

Win Ratio interim-timing analysis requires the primary validation caches
generated outside this repository. Pass their root explicitly:

```text
Rscript workflows/sensitivity/interim_timing/run_win_ratio_all_cases.R --output-dir /path/to/new_interim_results --formal-results-root=/path/to/primary_validation
```

For Win Ratio dependence analysis, prepare a compiled external runtime before
running the selected case and dependence level:

```text
Rscript workflows/sensitivity/dependence/prepare_win_ratio_runtime.R --output-dir /path/to/new_runtime
Rscript workflows/sensitivity/dependence/run_win_ratio_all_cases.R --output-dir /path/to/new_dependence_results --runtime-dir=/path/to/new_runtime --case=OAK --level=low --mode=formal
```

The Cox survival_colon calibration entry generates the complete H0 and
H1-standard PP paths for 10,000 attempts per scenario. It uses data seed 777
per scenario and model seed `12345 + replicate_id`. Exact-count screening and
paired-MCSE parameter selection are separate workflows. Independent validation
uses the separate streams recorded in `config/validation_seed_registry.csv`.
