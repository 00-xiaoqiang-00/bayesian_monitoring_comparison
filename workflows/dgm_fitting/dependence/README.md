# Adopted low/high-rho DGM refitting source

These workflows prepare observable PFS, OS and recurrence-CIF targets and fit the progression components at fixed dependence values. The frozen sensitivity parameters used by the monitoring analyses are in `config/rho_parameters.csv`.

The archived E07/E09 starting/fixed component parameters are numerically identical to the corresponding H1-standard/H1-observed rows in `config/dgm_parameters.csv`; this shared table is reused. Death margins and mixture weights stay fixed while progression parameters are fitted. H0 and H1-standard are derived using the locked standard-effect contrasts. The baseline-rho refit is a fit reference, not a replacement for the primary locked DGM.

From the code-package root:

```text
python -B workflows/dgm_fitting/dependence/prepare_km_targets.py --output-dir /path/to/new_external_results
python -B workflows/dgm_fitting/dependence/fit_rho_sensitivity_km_target.py --target-dir ../rho_targets --check-inputs
python -B workflows/dgm_fitting/dependence/fit_rho_sensitivity_km_target.py --target-dir ../rho_targets --output-dir /path/to/new_external_results
```

Use `--case OAK` (or another case) to limit a refit. All new output directories must be outside the source package. Dependencies are NumPy and SciPy, listed in the parent requirements file. The source performs deterministic fitting only; it does not run MCMC, monitoring simulations or threshold selection. Newly refitted results must be compared with the archived actual sensitivity table before use; SciPy/platform differences can affect optimizer rounding. The fitter reports its original absolute and baseline-degradation acceptance checks. The S1.3 manuscript curve-error grids are evaluated separately by `workflows/model_fit_assessment/`.
