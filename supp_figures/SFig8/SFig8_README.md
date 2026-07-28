# Supplementary Figure 8 — local reference-guided validation

This analysis addresses reviewer comment R1-2 using the existing Figure 3
Xenium breast-tumor tile and the archived SimSpace fitted parameter file. It
adds no dataset and modifies no package function.

## Scope

- Reference: one 1 mm × 1 mm Xenium tile (2,032 cells; 220 genes).
- Spatial evaluation: ten independently seeded SimSpace simulations generated
  from the same archived fitted parameter set (seeds 0–9).
- Molecular evaluation: the frozen seed-0 molecular realization aligned to
  the exactly replayed seed-0 spatial realization.
- Niche evaluation: the frozen BANKSY outputs already used for Figure 3.
- Inference: a local, within-source case study. Seeds, cells, cell types, and
  genes are not independent biological replicates, and no analysis is
  held-out.
- Calibration provenance: the archived Figure 3 parameter set was fitted with
  optimizer seed 0, population size 100, 30 generations, and one simulation
  replicate per fitness evaluation. The exact fitted file is
  `main_figures/Fig3/Panel_B_C_D_data/simspace_fitted_params.json`.

Moran's I is a calibration metric because the Figure 3 fit used cell-type
Moran's I and neighborhood entropy. Cell-type Geary's C, whole-layout
Ripley's L, gene-pattern metrics, permutation-centered cell-type
co-localization, and BANKSY niche summaries are out-of-objective diagnostics.

Cell-type Moran's I and Geary's C use a 20-nearest-neighbor graph and the seven
cell types represented by at least 20 cells in the reference tile. Gene
Moran's I and Geary's C use the 202 genes detected in 5%–95% of reference
cells. Ripley's L is intentionally not stratified by cell type or gene: it is
computed from all cell coordinates in the complete tile to evaluate overall
point-pattern scattering, matching Figure 2. Coordinates are scaled
independently along each axis to the unit square, and the uncorrected centered
profile L(d)−d is evaluated at 25 evenly spaced radii from 0 to 0.25.
Reference-versus-simulated agreement is computed across those 25 radii using
Pearson correlation and RMSE.

Spatial cell-cell communication inference is not included here. That distinct
question maps to R1-7 and requires a validated ligand–receptor truth model;
co-localization in this analysis is noncausal.

## Run

The frozen run used SimSpace 0.3.2 from commit
`ecf2855612871498dc89b8d43169229dfb8f6057`. From the reproducibility
repository root, activate the environment defined by `environment.yml`,
verify the adjacent source checkout, and run:

```bash
git -C ../SimSpace rev-parse HEAD
python -m pip install -e ../SimSpace
python supp_figures/SFig8/SFig8.py
```

The first command must report
`ecf2855612871498dc89b8d43169229dfb8f6057`.

The script writes:

- `supp_figures/SFig8/SFig8.png`;
- `example_output/SFig8/SFig8.png`; and
- frozen result and provenance tables under
  `supp_figures/SFig8/Panel_A_I_data/`.

`Panel_A_I_data/summary_metrics.tsv` is the numeric source for the manuscript
and response. `whole_layout_ripley_profiles.tsv` and
`whole_layout_ripley_agreement.tsv` contain the all-cell profiles and
ten-seed agreement values. `analysis_config.json`, `input_manifest.tsv`,
`seed_provenance.tsv`, and `software_versions.tsv` record the complete
configuration, input hashes, coordinate hashes, and execution environment.
