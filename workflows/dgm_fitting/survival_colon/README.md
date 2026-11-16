# survival_colon baseline DGM source

 entries retain the adopted likelihood/curve objectives and numerical settings. They use the paired recurrence/death records derived from the supplied public survival-colon patient table.

- `fit_joint_lognormal.R`: original paired-data censored bivariate-lognormal likelihood fit; includes dependence estimation and diagnostics.
- `fit_observed_curves.R`: the adopted observed-curve calibration, fitting OS first and recurrence cumulative incidence second while retaining the likelihood-based rho. The original rounded likelihood starts are retained. It outputs both the full-precision fit and the generator parameters using the original printed median/covariance precision.
- `derive_standard_h1.R`: original fixed-seed population-effect construction for H1-standard with a common control covariance. This performs 300,000 population draws per arm, not trial-replicate MCMC.

From the code-package root:

```text
Rscript workflows/dgm_fitting/survival_colon/fit_joint_lognormal.R --output-dir /path/to/new_external_results
Rscript workflows/dgm_fitting/survival_colon/fit_observed_curves.R --output-dir /path/to/new_external_results
Rscript workflows/dgm_fitting/survival_colon/derive_standard_h1.R --output-dir /path/to/new_external_results
```

R dependencies are mvtnorm, numDeriv, survival and MASS as used by the selected entry. Outputs must be outside the package. The observed-curve source uses the original rounded starting values; the first likelihood entry documents their origin. The historic observed-fit output, after the source's published rounding rule, matches the current frozen H1-observed generator. The formal simulation continues to read its locked configuration; a fresh optimizer result is not substituted automatically.
