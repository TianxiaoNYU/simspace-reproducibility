# Supplementary Figure 9 — local reference-guided validation

This analysis addresses reviewer comments R1-2, R1-3, R2-M6, and R3-2 using the
existing Figure 4 Xenium breast-tumor tile and the archived SimSpace fitted
parameter file. It adds no dataset and modifies no package function.

## Scope

- Reference: one 1 mm × 1 mm Xenium tile (2,032 cells; 220 genes).
- Spatial evaluation: ten independently seeded SimSpace simulations generated
  from the same archived fitted parameter set (seeds 0–9).
- Molecular evaluation: the frozen seed-0 molecular realization aligned to
  the exactly replayed seed-0 spatial realization, plus ten newly generated
  molecular realizations for spatial seeds 0–9 from one scDesign3 fit.
- Niche evaluation: the frozen BANKSY outputs already used for Figure 4.
- The reference and seed-0 cell maps are not repeated because they are already
  shown in Figure 4.
- Inference: a local, within-source case study. Seeds, cells, cell types, and
  genes are not independent biological replicates, and no analysis is
  held-out.
- Calibration provenance: the archived Figure 4 parameter set was fitted with
  optimizer seed 0, population size 100, 30 generations, and one simulation
  replicate per fitness evaluation. The exact fitted file is
  `main_figures/Fig4/Panel_B_C_D_data/simspace_fitted_params.json`.

Moran's I is a calibration metric because the Figure 4 fit used cell-type
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
Pearson correlation and RMSE. Panel A compares every simulated seed directly
with the Xenium reference; it is not a simulation-versus-simulation
correlation.

Expression fidelity uses all 220 genes shared by the Xenium reference and
simulated datasets. Panels H and I compare the frozen seed-0 realization with
Xenium using the log1p gene-wise raw-count mean and log1p unbiased raw-count
variance. Panels J and K fit the scDesign3 marginal and copula models once,
generate new molecular profiles for spatial seeds 0–9 using the matching
molecular seeds, and separately report PCC and RMSE between each realization's
log1p gene-mean vector and the Xenium vector across all 220 genes. These
same-tile results measure reference concordance rather than held-out accuracy.

Spatial cell-cell communication inference is not included here. That distinct
question maps to R1-7 and requires a validated ligand–receptor truth model;
co-localization in this analysis is noncausal.

## Run

The frozen run used SimSpace 0.3.4 from commit
`de0a4c002e4ae733e354e3e180ab69b381ad994a`. From the reproducibility
repository root, activate the environment defined by `environment.yml`,
verify the adjacent source checkout, and run:

```bash
git -C ../SimSpace rev-parse HEAD
python -m pip install -e ../SimSpace
python supp_figures/SFig9/SFig9.py
```

The first command must report
`de0a4c002e4ae733e354e3e180ab69b381ad994a`.

To regenerate Panels J–K, use R 4.4.2 with scDesign3 1.5.0 and
SingleCellExperiment 1.28.1:

```bash
Rscript supp_figures/SFig9/Panel_L_src/generate_molecular_replicates.R \
  main_figures/Fig4/Panel_B_C_D_data/Xenium_reference_metadata.csv \
  main_figures/Fig4/Panel_B_C_D_data/Xenium_reference_count.csv \
  supp_figures/SFig9/Panel_L_data/molecular_simulation_design.tsv \
  supp_figures/SFig9/Panel_L_data/molecular_replicate_summaries.tsv
python supp_figures/SFig9/SFig9.py
```

The initial Python run can use the frozen molecular summary supplied in
`Panel_L_data/`; the R command is needed only to regenerate that table.

The script writes:

- `supp_figures/SFig9/SFig9.png`;
- `example_output/SFig9/SFig9.png`; and
- frozen result and provenance tables under
  `supp_figures/SFig9/Panel_A_I_data/`;
- seed-0 expression-fidelity tables under
  `supp_figures/SFig9/Panel_J_K_data/`; and
- ten-seed molecular design, summaries, and agreement under
  `supp_figures/SFig9/Panel_L_data/`.

`Panel_A_I_data/summary_metrics.tsv` is the numeric source for the manuscript
and response. `whole_layout_ripley_profiles.tsv` and
`whole_layout_ripley_agreement.tsv` contain the all-cell profiles and
ten-seed agreement values. `Panel_J_K_data/expression_fidelity.tsv` and
`expression_agreement.tsv` contain the 220-gene seed-0 results.
`Panel_L_data/molecular_replicate_summaries.tsv` and
`molecular_replicate_agreement.tsv` contain the ten-seed molecular results,
and `molecular_software_versions.tsv` records the R environment.
`analysis_config.json`, `input_manifest.tsv`, `seed_provenance.tsv`, and
`software_versions.tsv` record the configuration, input hashes, coordinate
hashes, and Python execution environment; the R versions for Panels J–K are
specified above.
