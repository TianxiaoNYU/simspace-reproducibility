# Supplementary Figure 7 — reference-free 3D simulation

This analysis addresses reviewer comment R1-5 with a focused quantitative
evaluation of SimSpace's current three-dimensional support. It evaluates one
fixed reference-free configuration and does not use or claim a
reference-guided 3D model.

## Design

- Ten simulation seeds: 0–9.
- Volume: 50 × 50 × 20 lattice sites.
- Spatial model: three MRF niches, six phenotypes, a six-connected 3D
  neighborhood, spatial smoothness parameter 5, ten niche sweeps, and four
  phenotype sweeps.
- Cell output: 0.85 uniform density retention and coordinate perturbation
  standard deviation 0.40 lattice units, yielding approximately 42,500
  retained cells per seed.
- Molecular output: 36 Gamma–Poisson genes, including 24 phenotype marker
  genes and 12 background genes. Counts are conditioned on phenotype; no
  direct coordinate-conditioned expression effect is used.

The analysis asks only two quantitative questions:

1. Do the simulated niches, phenotypes, and marker-gene profiles exhibit
   spatial organization in the complete volume? Three-dimensional Moran's I
   is computed on an 18-nearest-neighbor graph. Within each seed, niche and
   phenotype results are averaged across indicators, while marker-gene
   results are summarized by the median across marker genes.
2. Can the exported niche truth be recovered in a simple downstream
   analysis? K-means clustering is applied to either standardized
   single-cell log-counts alone or those features concatenated with a
   separately standardized mean log-count vector for each cell and its 24
   nearest neighbors in 3D. Recovery is evaluated with ARI and NMI against
   the exported niche label.

Seeds and within-seed feature summaries are simulation replicates and
components, not independent biological replicates. The expression-only
clustering is a baseline for the stated downstream task, not a comparison
with another simulator or a published domain-identification method.

## Fast figure regeneration

The archived run uses SimSpace 0.3.4 from commit
`de0a4c002e4ae733e354e3e180ab69b381ad994a`. From the reproducibility
repository root, the normal entry point reads the archived ten-seed outputs
and regenerates only the figure:

```bash
conda activate simspace-repro
python supp_figures/SFig7/SFig7.py
```

This fast path writes:

- `supp_figures/SFig7/SFig7.png`;
- `example_output/SFig7/SFig7.png`.

It does not rerun SimSpace or overwrite the archived metrics, cell metadata,
or molecular counts. Rendering takes approximately 5 seconds after the
Matplotlib font cache is initialized; a clean environment may spend an
additional 10--20 seconds building that cache on the first run.

## Full simulation and data regeneration

To rerun all ten simulations and refresh every archived result:

```bash
conda activate simspace-repro
python supp_figures/SFig7/SFig7_generate_data.py
```

The full generation script takes approximately 6--7 minutes in the
development environment and writes:

- both figure copies listed above;
- `Panel_A_F_data/spatial_structure_metrics.tsv`, component-level Moran's I;
- `Panel_A_F_data/spatial_structure_summary.tsv`, the seed-level spatial
  summary used in the figure;
- `Panel_A_F_data/domain_recovery_metrics.tsv`, the frozen ARI/NMI results;
- `Panel_A_F_data/summary_metrics.tsv`, the single seed-level numeric source
  for all values reported in the manuscript and response;
- compressed cell metadata and molecular counts for all ten seeds;
- the complete analysis configuration and software versions.

## Interpretation

The figure evaluates volumetric organization and one truth-based downstream
use under the declared reference-free configuration. It does not establish
agreement with a real 3D tissue, support reference-guided 3D calibration, or
imply a direct coordinate-conditioned molecular effect.

Across the ten seeds, the median three-dimensional Moran's I was 0.592 for
niche indicators (range 0.574–0.611), 0.247 for phenotype indicators
(0.241–0.254), and 0.165 for marker-gene expression (0.157–0.189).
Adding three-dimensional neighborhood features improved niche recovery in
every seed: median ARI increased from 0.219 (0.142–0.591) to 0.811
(0.793–0.833), and median NMI increased from 0.363 (0.254–0.649) to 0.735
(0.717–0.761). Across all ten seeds, the median paired gains were 0.594 for
ARI and 0.364 for NMI.
