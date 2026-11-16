# IMpower131 adopted baseline DGM fit

Run from the package root:

```sh
Rscript workflows/dgm_fitting/impower131/fit_impower131_baseline.R --output-dir /absolute/external/impower131
```

An explicit output directory outside the source package is required. A single
positional directory is also accepted. Relative paths resolve from the current
working directory; Windows paths such as `D:/analysis_runs/impower131` work.
Missing output arguments and destinations within the source package are rejected
before any fitting or output. R and the `survival` package are required.

The entry fits the death margins first, retains their selected parameters,
then fits latent progression at the adopted rho = 0.50. It preserves the
original candidate models, minimax/MLE selection, median constraints, common
random numbers, objective, optimization steps and diagnostics. This entry
performs numerical calibration; it does not read a finished
parameter table as a fitted result.

The fixed input is `inputs/reconstructed_data/fitting_anchors/impower131/IMpower131_marginal_ipd.rds`.
It contains only the four historical marginal pseudo-IPD fields actually read
by this fit (`ipd$os_trt$IPD`, `os_ctrl`, `pfs_trt`, `pfs_ctrl`). Time is in
months; status 1 is an event and 0 is censoring. The original values are
preserved. They are the earlier baseline-fitting anchors and differ from
the separately reconstructed 2026-07-31 validation pseudo-IPD. PFS and OS
records must not be paired by row.

The historical rho = 0.50 output, after its original four-decimal export
of exp(mu) and covariance, gives control exp(mu) = (6.1202, 13.4685),
Sigma = ((0.4209, 0.3287), (0.3287, 1.0268)); observed treatment exp(mu) =
(7.5610, 14.3431), Sigma = ((0.8207, 0.5298), (0.5298, 1.3679)). These are
the adopted E05-v6 baseline parameters. The standard alternative is a
subsequent design construction with exp(mu) = (8, 18) and the control
covariance; it is not an observed-curve refit.

The original fit's tolerance flags are retained. Its report marked no rho
as satisfying every earlier main-simulation criterion. The later E05 curve
assessment uses the adopted rounded model over PFS 0–28 and OS 0–36 months;
those evaluation horizons do not change this historical optimization.
