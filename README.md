# Bayesian monitoring comparison: GitHub code package

This repository contains the source code, locked configurations and fixed public
inputs needed to reproduce the Bayesian Win Ratio, Cox, restricted mean survival
time (RMST) and restricted mean utility time (RMUT) comparison in the IMpower131,
KEYNOTE-598, OAK and survival_colon settings.


| Directory | Purpose |
|---|---|
| `config/` | Locked DGM parameters, trial processes, thresholds, seeds and calculation settings |
| `inputs/` | Public clinical data, reconstructed marginal data and fixed numerical inputs |
| `src/` | Reusable data-generating, statistical and monitoring implementations |
| `workflows/` | Analysis entry points grouped by task |
| `third_party/` | Included Win Ratio package sources and license notices |
| `docs/` | Runtime instructions and third-party attribution |

Start with `analysis_index.csv` to find the workflow for a method, case and
analysis. Read `docs/RUNNING.md` before running a formal workflow. All generated
caches, fitted objects, posterior draws, compiled libraries and run outputs must
be written outside this repository.

RMUT parameter selection shows help when run without arguments. Use `rescreen`
to regenerate candidate and selected tables, or `verify --reference-dir DIR`
to check externally supplied tables. See
[`workflows/parameter_selection/rmut/README.md`](workflows/parameter_selection/rmut/README.md).

The eight calibration PP-path RDS files and two historical marginal-IPD fitting
anchors under `inputs/` are fixed analysis inputs. They are deliberately retained
and listed in the input manifests. Reconstructed PFS and OS records are marginal
inputs and are not matched patient pairs. The survival_colon input is derived from
the public `survival` package.

The population-estimand verification workflows can verify an external results
archive when one is supplied; result tables are not bundled here. Explicit
recomputation writes to a new external output directory. Figure-generation
workflows and figure-only inputs are omitted from this GitHub package by design.

`FILES_SHA256.csv` records the distributed files except the manifest itself.
The third-party attribution and permission scope are documented in
`docs/THIRD_PARTY_NOTICES.md`.
