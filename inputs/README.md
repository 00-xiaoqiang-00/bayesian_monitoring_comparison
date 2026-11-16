# Fixed analysis inputs

| Directory | Contents |
|---|---|
| `clinical_data/` | Digitized PFS/OS coordinates, numbers at risk, source definitions and public survival_colon records |
| `reconstructed_data/marginals/` | Separate reconstructed PFS and OS marginal samples |
| `reconstructed_data/fitting_anchors/` | Marginal-IPD inputs used by the adopted IMpower131 and OAK baseline fits |
| `calibration_pp/` | Eight fixed H0/H1-standard PP-path inputs for RMUT selection, with an input manifest |

The baseline-fitting anchors and validation marginal reconstructions are
different numerical inputs. Select the inputs required by the relevant
workflow. Separate PFS/OS marginal samples do not identify patient pairs.

`clinical_data/source_manifest.csv` describes the clinical sources and
coordinate identities. It also records the OAK month-3 OS risk-set correction
from the FDA source. The accompanying risk-set and evaluation tables define
the supplied targets and evaluation grids.

The survival_colon records are derived from the public `survival::colon`
dataset. The included clinical-trial curves and reconstructed marginals are
not restricted individual clinical-trial records.

Figure-only inputs and plotting assets are intentionally omitted from this
GitHub code package. The full submission distribution contains the fixed Figure
S1 reference tables when those materials are needed.
