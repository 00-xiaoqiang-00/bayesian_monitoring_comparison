# Cox implementation

This directory contains the statistical implementation and case-specific
drivers for Cox. Run the public entries under `workflows/`, which provide
the fixed configuration, source paths and an external output location.

`--preflight-only` checks input/source availability and R syntax without model
loading, simulation or MCMC. A numerical run requires the R dependencies and,
for Stan or third-party C++ routines, a working compiler toolchain.

Locked thresholds, trial settings and random streams are supplied in `config/`.
Model outputs, fitted objects and PP caches must be written outside the whole
repository. Refer to the root entry-point index for calibration,
independent validation and sensitivity commands.

## survival_colon calibration PP inputs

`workflows/calibration/run_cox_survival_colon.R` generates the two full
H0/H1-standard PP caches for this case. Each scenario has 10,000 attempted
records, including unsuccessful fits. Data generation uses seed 777 per
scenario; posterior fitting uses seed 12345+j. The model, diagnostic filter,
trial process and precompilation model are defined by the calibration source
identified in `config/cox_colon_calibration.json`.

This entry stops after saving the PP paths. Parameter eligibility and final
selection are separate workflows using exact H0 500/10,000 success and the
paired-MCSE rule. No early-stopping replay from the source script is performed
by the PP-generation entry. H1-observed does not enter calibration selection.
