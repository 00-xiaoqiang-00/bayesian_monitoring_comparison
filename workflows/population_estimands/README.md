# Population estimands

The three entries evaluate the locked DGM under H1-standard, H1-observed and
the scheduled interim/final restriction windows. Calculation and verification
settings are under `config/estimands/`. Archived result tables are supplied
separately with the full submission distribution.

From the code-package root:

```text
python -B workflows/population_estimands/compute_dgm_estimands.py --verify-archive --archive /path/to/external_archive/dgm_estimands.csv
python -B workflows/population_estimands/compute_h1_observed_estimands.py --verify-archive --archive-dir /path/to/external_archive
python -B workflows/population_estimands/compute_interim_estimands.py --verify-archive --reference-csv /path/to/external_archive/interim_true_estimands.csv
```

These verification modes check saved numeric results and locked inputs without
model fitting, patient generation, or new Monte Carlo integration. They require
an external results archive; a standalone code directory may locate archives
explicitly with `--archive`, `--archive-dir`, or `--reference-csv`, respectively.

Explicit `--recompute --output-dir /path/to/new_external_directory` requests a
calculation. H1-observed additionally supports `--integrals-only`, which only
evaluates deterministic RMST/RMUT and does not replace the full estimand table.
H1-standard uses 12 Cox batches of 100,000 observations per arm and 20 WR
batches of 200,000 pairs. H1-observed uses 24 Cox batches of 250,000 observations
per arm and 100 WR batches of 200,000 pairs, preserving the configured streams.
The interim calculation uses 1,500 trapezoidal intervals (1,501 nodes).

Cox targets the marginal first-event projection under the locked final trial
process. WR targets the death-prioritized IPCW weighted win/loss ratio. RMST and
RMUT are contrasts in months. KEYNOTE-598 retains separate final windows:
12 months for RMST and 18 months for RMUT. Numerical uncertainty excludes
original-trial sampling and DGM-fit uncertainty. Scientific configurations
record the relevant numerical environments and fixed computation settings.
