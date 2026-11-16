# Frozen-DGM curve diagnostics

This entry evaluates the formal primary and actual low/high-rho
H1-observed data-generating models for OAK, IMpower131, KEYNOTE-598 and
survival_colon. It calculates MAE, RMSE, risk-set-weighted RMSE and maximum
absolute curve error from frozen parameters. It does not fit parameters,
run monitoring replicates, or change the package inputs.

## Run

Use Python 3.10 or newer, NumPy and SciPy:

```text
python -m pip install -r requirements.txt
python -B workflows/model_fit_assessment/evaluate_rho_curve_fit.py --output-dir /path/to/new_external_results
python -B workflows/model_fit_assessment/evaluate_rho_curve_fit.py --case KEYNOTE-598 --output-dir /path/to/new_external_results --write-curves
```

The package root is located from the script, so the current working
directory need not be the repository root. `--package-root` selects another
complete copy of the package. Choose an output directory outside the
package. All generated diagnostics are written there.

## Parameters and targets

- Primary models: `config/dgm_parameters.csv`, selecting
  H1-observed. Effective correlations from the frozen covariance matrices
  are retained. The primary colon correlations differ between arms.
- Actual sensitivity models:
  `config/rho_parameters.csv`, selecting
  H1-observed and low/high. Its baseline rows are not substituted for the
  formal primary models. Mixture weights and death margins remain at the
  values in the actual table.
- Lung-cancer targets: current digitized coordinates and independent
  marginal pseudo-IPD in `inputs/reconstructed_data/marginals/`. The marginal records
  are used for risk counts, not patient-level endpoint pairing.
- Colon target: `inputs/clinical_data/survival_colon.csv`.
  Obs and Lev+5FU comprise 619 patients, with recurrence and death records
  joined by patient id. Days are converted with 30.4375 days/month.

The target curve and risk counts are reconstructed directly from those
inputs. No local historical directory or archived source script is needed.

## Evaluation grids and curve definitions

Each lung-cancer curve has 1,000 evaluation points from 0.03 months to its
right endpoint minus 0.0001 months, matching the S1.3 numerical grid.
Digitized targets use linear interpolation. The displayed windows start
at zero, while this original numerical grid starts at 0.03 months.

| Case | PFS right endpoint | OS right endpoint |
|---|---|---|
| OAK | 18 months | 27 months |
| IMpower131 | 28 months | 36 months |
| KEYNOTE-598 | Per endpoint and arm | Per endpoint and arm |

For KEYNOTE-598, each endpoint-by-arm curve has its own right endpoint:
the earlier of the last digitized coordinate and the fifth-largest
follow-up time in its marginal pseudo-IPD. The four endpoints are
26.730000/25.638750 months for control/treatment PFS and
29.393077/28.611000 months for control/treatment OS. Every evaluation point
has at least five subjects at risk. All four error measures use the same
grid for that curve, reused across the primary, low-rho and high-rho models.

Lung-cancer PFS is Pr(T_R > t, T_D > t), evaluated by Gaussian quadrature;
OS is Pr(T_D > t). The default is 640 quadrature nodes. No reoptimization
is performed at another rho.

Colon uses 169 points at 0.5-month intervals over 0–84 months. The observed
recurrence curve is an Aalen–Johansen CIF with death as a competing event;
same-day death takes precedence. Observed OS uses the death-record KM
curve. The primary and actual low/high-rho models all use analytic curves:
recurrence CIF is Pr(T_R <= t, T_R < T_D), and OS is Pr(T_D > t). The
generator does not enforce T_R <= T_D.

For every curve, weighted RMSE uses
`max(nrisk(t)/nrisk(start), 0.025)`. Ordinary MAE and RMSE use equal weights
over the evaluation grid. Models within one case, arm and endpoint use
identical targets, grids and risk counts.

## Fixed references and outputs

The finite-trial Figure S1 reference table and the archived 48-row metrics
table are submission materials and are intentionally omitted from this GitHub
code package. The diagnostic still evaluates the analytic primary and
low/high-rho curves from the fixed clinical inputs. For the finite-trial
Monte Carlo comparison, use the corresponding workflow and references from the
full submission distribution. This GitHub entry does not load figure-only
inputs.

The default output includes one metrics table, one evaluation-window table,
the parameters actually used, and primary-reference comparisons for the
lung-cancer curves using `inputs/clinical_data/evaluation_curve_sources.csv`.
Colon assessment uses the analytic curves described above. `--write-curves`
also saves the pointwise target, model and risk counts. The required output
directory is external; generated figures, diagnostic logs and fitted-model
objects are not stored in this source package.
