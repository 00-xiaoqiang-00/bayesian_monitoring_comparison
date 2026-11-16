# DGM fitting and construction

These entries provide the adopted fitting and construction steps used before the monitoring analyses. They retain the original numerical objectives and selection rules. Formal monitoring continues to read the frozen configuration in `config/`; a fresh fit never overwrites it automatically.

| Directory | Adopted workflow |
|---|---|
| `impower131/` | Original rho=0.50 OS-first baseline calibration using its historical marginal-IPD anchors; original export precision recovers the adopted E05-v6 parameters |
| `oak_os_margin/` | Original OS death-margin candidate fitting and selection from the historical OS-IPD anchors |
| `oak/` | Formal two-component progression fitting with the adopted death margins locked; selected local least-squares candidate is retained |
| `keynote598/` | Shared-component mixture weight and PFS/OS parameters jointly fitted, with within-component rho=0.50 |
| `survival_colon/` | Original paired-data likelihood, observed OS/CIF calibration, and H1-standard construction |
| `dependence/` | Actual low/high-rho observable-target refitting; death margins and mixture weights stay fixed, followed by the locked standard-effect construction |

Read the case README before execution. The IMpower131/OAK baseline-fitting anchor IPD differ from the marginal validation reconstructions and are supplied separately under `inputs/reconstructed_data/fitting_anchors/`. They contain raw marginal inputs only, not fitted model objects or machine configuration. Separate PFS/OS marginal records must not be paired row by row.

The Python entries require Python 3.10+ and NumPy/SciPy from `requirements.txt`. IMpower131/OAK OS entries require R and survival; the selected colon entry also uses mvtnorm, numDeriv or MASS as listed in its source. All outputs must go outside the code package.

Examples from the code-package root:

```text
Rscript workflows/dgm_fitting/impower131/fit_impower131_baseline.R --output-dir /path/to/new_external_results
Rscript workflows/dgm_fitting/oak_os_margin/fit_oak_os_margin.R --output-dir /path/to/new_external_results
python -B workflows/dgm_fitting/oak/fit_oak_dgm.py --output-dir /path/to/new_external_results
python -B workflows/dgm_fitting/keynote598/fit_keynote598.py --output-dir /path/to/new_external_results --require-reference-match
Rscript workflows/dgm_fitting/survival_colon/fit_observed_curves.R --output-dir /path/to/new_external_results
python -B workflows/dgm_fitting/dependence/prepare_km_targets.py --output-dir /path/to/new_external_results
python -B workflows/dgm_fitting/dependence/fit_rho_sensitivity_km_target.py --target-dir ../rho_targets --output-dir /path/to/new_external_results
```

Baseline and dependence fitting are distinct from decision-boundary calibration, candidate screening and independent validation. Low/high-rho curve-error evaluation on the manuscript S1.3 grids is in `workflows/model_fit_assessment/`. The original dependence-refit acceptance checks are retained in the fitting source and should not be confused with the later diagnostic evaluation grids.
