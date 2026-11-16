# Marginal pseudo-IPD reconstruction

These entries reproduce the adopted marginal reconstructions using the supplied digitized curves and the published risk/event totals embedded in the source.

- `reconstruct_oak_impower_marginals.py`: OAK and IMpower131, using the marginal reconstruction algorithm.
- `reconstruct_keynote598_marginals.py`: KEYNOTE-598, using only the marginal-reconstruction portion of the original reconstruction-and-fit source. 

From the code-package root:

```text
python -B workflows/reconstruction/reconstruct_oak_impower_marginals.py --output-dir /path/to/new_external_results
python -B workflows/reconstruction/reconstruct_keynote598_marginals.py --output-dir /path/to/new_external_results
```

Choose new output directories outside the package. Inputs default to `inputs/clinical_data/`. No MCMC or trial-replicate simulation is run. PFS and OS are separate marginal reconstructions and must not be paired row by row.

The frozen OAK OS treatment risk count at month 3 is 363, verified against Figure 1 of the FDA CDTL review (PDF page 13; printed page 12). Lancet Figure 2A prints 263 at this position. The FDA source URL and location are recorded in `inputs/clinical_data/source_manifest.csv`. The final tail interval accounts for subjects remaining beyond the last printed risk count.
