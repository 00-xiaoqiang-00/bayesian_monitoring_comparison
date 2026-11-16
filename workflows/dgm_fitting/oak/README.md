# OAK: fit the adopted progression mixtures with OS margins fixed

This entry point reproduces the deterministic curve-fitting step that produced
the adopted OAK `H1-observed` progression-mixture parameters at `rho=0.70`.
It performs numerical optimization; the locked reference parameters are read
only after the fit for comparison.

## Model and order of fitting

Within each arm, the latent progression time has a two-component Lognormal
mixture. Both components share the same arm-specific death-time distribution.
PFS is derived as `min(TR, TD)`.

1. Read the existing OS death-margin lock. The control values are
   `mu_D=2.3088355183235585, sigma_D=1.0177`; treatment uses
   `mu_D=2.636754680552583, sigma_D=1.1709`. These values came from the
   earlier IPD-anchored OS fitting record and are not re-estimated here.
2. Fit the control progression mixture to the control PFS curve, its median
   anchor, and the 6-, 12-, and 18-month survival anchors. The adopted control
   fit used a progression-component standard-deviation lower bound of `0.30`.
3. Keep the fitted control distribution and both death margins fixed. Fit
   the treatment progression mixture to the treatment PFS curve, median,
   landmark survival probabilities, and a curve-level PFS Cox HR target.
4. Report the fitted parameters and their differences from the archived
   formal parameters. Do not overwrite the study's locked DGM configuration.

Only five variables are optimized per arm:

```text
pi, mu_R1, log(sigma_R1), log(mu_R2-mu_R1), log(sigma_R2)
```

There are no free death parameters or OS curve residuals in these objectives.
The fixed correlation is a modeling assumption; it is not identified from
the separate marginal curves.

## Candidate adoption

The adopted treatment result was the **first local least-squares candidate**
started from the earlier treatment-mixture fit. Its archived weight is
`pi=0.2623242204350347`. The corresponding control weight is
`pi=0.13123216948698738`.

A subsequent historical SLSQP refinement produced another candidate with
`pi=0.26734832311797696`, but that candidate was not written to the formal
configuration. It is not used by this entry point. The selection is preserved
as an explicit historical adoption decision, not replaced by choosing the
smallest numerical difference among near-equal optimizer results.

The initial treatment point is the earlier fit with `pi=0.2264164067590971`;
it is not the answer taken from the reference file. The legacy standalone
`oak_mixture_calibration.py` described the earlier single-control/mixture-
treatment version and is not a complete entry point for the adopted two-arm
mixture model. The relevant numerical kernel has been retained in
`oak_fit_core.py`, including its 40-node quadrature and curve-level HR
calculation, to preserve the original objective.

## Run

Python 3.10 or later, NumPy, and SciPy are required. Use the shared DGM fitting
requirements file. From the package root, for example:

```bash
python workflows/dgm_fitting/oak/fit_oak_dgm.py --output-dir /path/to/new_external_results
```

`--output-dir` is required and must be outside the package and the input
directories. The command neither updates locked parameters nor launches any
trial simulation or MCMC.

Defaults resolve relative to the package location, so the command also works
from a different current directory. Explicit overrides are available:

```bash
python workflows/dgm_fitting/oak/fit_oak_dgm.py \
  --curve-dir inputs/clinical_data/OAK \
  --input-dir config/dgm_fitting/oak \
  --output-dir /path/to/new_external_results
```

The curve directory supplies the existing headerless two-column PFS CSVs.
The metadata directory supplies `oak_locked_death_parameters.json`,
`fit_settings.json`, and `locked_reference.csv`.

## Results and numerical comparison

The command writes:

- `fitted_parameters.csv`: newly optimized arm-specific parameters.
- `parameter_comparison.csv`: fitted versus archived parameters.
- `fitted_reference_curves.csv`: PFS and OS curves on the comparison grid.
- `fit_report.json`: optimizer convergence, curve diagnostics, candidate
  adoption, software versions, tolerances, and the comparison result.

The default absolute tolerances are `5e-5` for parameters and `1e-5` for
survival probabilities. These are numerical refitting checks, not changes to
the full-precision archived parameters. The command exits with code `2` if
the fitted-reference comparison or the historical PFS acceptance checks fail.
Such a failure does not replace any locked DGM input.

The reported PFS Cox HR is the original deterministic curve diagnostic through
24 months. It is not the paper's final-process Cox HR projection. The module
does not reselect `rho`, construct `H1-standard`, rerun dependence sensitivity,
or reproduce the earlier marginal pseudo-IPD reconstruction.
