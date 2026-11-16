# KEYNOTE-598 baseline DGM curve calibration

This entry point estimates the shared component weight and the progression
and death component parameters jointly against each arm's PFS and OS curves.
One component label selects both endpoint parameters. The within-component
log-time correlation is fixed at 0.5. Separate marginal curves do not identify
patient-level pairing, component membership or the dependence parameter.

## Run

From the repository root, with NumPy and SciPy installed:

```text
python workflows/dgm_fitting/keynote598/fit_keynote598.py --output-dir /absolute/path/outside/package/keynote598_fit
```

The default initialization/reference inputs are
`config/dgm_fitting/keynote598/` relative to the repository root.
The source curves reuse `inputs/clinical_data/KEYNOTE-598/`.
`--input-dir` and `--curve-dir` can select alternative input copies.
The output directory must be new and outside this package. The program never
changes `config/` or the locked parameters used by the simulations.
An output directory inside either custom input directory is also rejected.
`--arm treatment` or `--arm control` runs one arm. By default both are fitted.

For a numerical replication check, add `--require-reference-match`. The
default parameter comparison tolerance is 1e-4 and may be set explicitly with
`--parameter-tolerance`. This comparison is performed after estimation and
does not alter initialization, objective, constraints or candidate selection.

## Numerical specification

- Nine free coordinates encode pi, the two ordered progression locations and
  scales, and the two ordered death locations and scales.
- Pi is bounded to [0.05, 0.95]; progression log standard deviations are at
  least 0.15; death log standard deviations are at most 1.5.
- Each arm uses six initial values. Four combine parameters from earlier
  numerical candidate models; two are generic starts. The supplied earlier
  candidate tables are initialization inputs only, not fitted output from
  this command and not clinically identified patient classes.
- The objective concatenates 520-point weighted PFS and OS curve residuals,
  each divided by 0.02, and median/12-month survival anchors divided by 0.008
  and repeated eight times. Risk weights have the original 0.005 floor.
- The PFS probability is evaluated with the original 48-node Gauss-Legendre
  integration. Each initial value uses bounded least squares with 18,000
  maximum evaluations and xtol/ftol/gtol=1e-11.
- SLSQP refines the best local candidate (700 iterations, ftol=1e-12).
  Constraints require maximum supported-curve difference <=0.05, median
  difference <=0.25 months and 12-month survival difference <=0.012.
  The original feasibility tolerance is 2e-5. Candidate selection preserves
  the original feasibility/constraint-violation/objective ordering.


## Inputs and outputs

The four `KEYNOTE598_*` files contain the digitized PFS/OS curve coordinates
from the source trial's published curves. Risk tables used in weighting are
included explicitly in the script. The two `initialization_*` tables supply
the earlier candidate fits described above. `reference_parameters.csv` stores
the selected regularized fit for post-fit comparison only. It is never used
as an optimizer start.

The fitted parameters, optimizer diagnostics, parameter comparison, curve
comparison and runtime/input hashes are written as new CSV/JSON files.
Numerical optimizer results may differ slightly across NumPy/SciPy versions;
both parameter and survival-curve differences are reported. These fits do not
automatically replace the locked DGM.

The study uses rounded locked values pi_1=0.1936579096 for treatment and
pi_1=0.2070545742 for control; the reference table retains the original
full-precision numerical fit. Mixture components serve numerical curve
reproduction and are not clinical subgroups. This routine does not perform
MCMC, operating-characteristic simulations or rho sensitivity analysis.
