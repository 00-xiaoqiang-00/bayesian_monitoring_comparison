# RMUT implementation

This directory contains the statistical implementation and case-specific
drivers for RMUT. Run the public entries under `workflows/`, which provide
the fixed configuration, source paths and an external output location.

`--preflight-only` checks input/source availability and R syntax without model
loading, simulation or MCMC. A numerical run requires the R dependencies and,
for Stan or third-party C++ routines, a working compiler toolchain.

Locked thresholds, trial settings and random streams are supplied in `config/`.
Model outputs, fitted objects and PP caches must be written outside the whole
repository. Refer to the root entry-point index for calibration,
independent validation and sensitivity commands.
