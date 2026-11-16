# OAK death-margin source fit

Run from the package root:

```sh
Rscript workflows/dgm_fitting/oak_os_margin/fit_oak_os_margin.R --output-dir /absolute/external/oak_os_margin
```

An explicit output directory outside the source package is required. A single
positional directory is also accepted. Relative paths resolve from the current
working directory; Windows paths such as `D:/analysis_runs/oak_os_margin` work.
Missing output arguments and destinations within the source package are rejected
before any fitting or output. R and the `survival` package are required.

This entry retains the original OS-first death-margin fit. It uses each
reconstructed OS KM median as the anchor, considers fixed-median MLE,
weighted KM, minimax KM, flexible-median minimax and clipped/unrestricted
MLE candidates, and selects them using the original curve-error rule.
The median window is ±2%, the main fitting window is the 90th percentile
of observed times, and the curve tolerance is 0.03.

The fixed input is `inputs/reconstructed_data/fitting_anchors/oak_os_margin/OAK_os_marginal_ipd.rds`.
It contains only `ipd$os_trt$IPD` and `ipd$os_ctrl$IPD`, preserving the
historical baseline-fitting values. Time is in months; status 1 is death
and 0 is censoring. These anchors differ from the 2026-07-31 validation
pseudo-IPD. The input contains no fitted parameters or machine configuration.

The adopted death margins use the recorded export precision:
control exp(mu_D) = 10.0627 and sigma_D = 1.0177; observed treatment
exp(mu_D) = 13.9678 and sigma_D = 1.1709. Both historical selected methods
were flexible-median minimax KM. These margins are shared by the two
progression components in the formal OAK DGM. The formal standard
alternative subsequently uses control sigma_D = 1.0177 in both arms.

Only the OS step is published here. The earlier single-component PFS
model in the source project is not the formal two-component OAK model.
The progression mixture fitting entry remains in the adjacent `oak` directory.
