# Supplementary Figure 9 — local reference-guided validation

This analysis addresses reviewer comments R1-2, R1-3, R1-8, R2-M6, and R3-2
using the Figure 4 Xenium breast-tumor tile. It modifies neither SimSpace
package functions nor Figure 4 or any other figure.

## Scope and design

- Reference: one 1 mm × 1 mm Xenium tile (2,032 cells; 220 genes).
- Independent spatial fitting: ten archived fitted-parameter JSONs obtained
  from fully independent genetic-algorithm calibrations with optimizer seeds
  0–9. The author-confirmed original-file mapping is suffix 1–10 to optimizer
  seeds 0–9. The files are stored with unambiguous zero-padded names under
  `Panel_A_I_data/fitted_params/`.
- Spatial generation: one layout is generated per fitted parameter set using
  the matching generation seed 0–9. State labels from each independent fit are
  mapped to Xenium cell types by abundance rank (`state_rank`), rather than by
  the permutation-arbitrary raw state identifier.
- Molecular robustness: scDesign3 is fitted once to the Xenium counts, then
  one molecular draw is generated for each independently fitted layout using
  the matching molecular seed 0–9.
- Frozen Figure 4 diagnostics: Panels B–C and G–I retain the original Figure 4
  spatial/count realization and BANKSY outputs as a separate, explicitly
  labeled diagnostic. They are not presented as results from the ten new
  calibrations, and Figure 4 itself is unchanged.
- Inference: this is a local, same-tile case study, not held-out validation.
  Runs, cells, cell types, and genes are not independent biological
  replicates.

`Panel_A_I_data/calibration_manifest.tsv` records the optimizer,
generation, and molecular seed mapping and the fitted-file paths, hashes, and
sizes. `Panel_A_I_data/seed_provenance.tsv` additionally records generated
cell counts and coordinate hashes. The script validates the ten JSON schemas,
their 83 finite fitted values, seeds 0–9, and ten unique parameter hashes. It
consumes the archived fitted outputs; it does not rerun the genetic-algorithm
calibrations.

Moran's I is a calibration metric because the spatial fit used cell-type
Moran's I and neighborhood entropy. Cell-type Geary's C, whole-layout
Ripley's L, gene-pattern metrics, permutation-centered cell-type
co-localization, and BANKSY niche summaries are out-of-objective diagnostics.
Cell-type metrics use a 20-nearest-neighbor graph and the seven cell types with
at least 20 reference cells. Gene spatial metrics use the 202 genes detected
in 5%–95% of reference cells. Ripley's L uses all coordinates, after separate
axis normalization to the unit square, at 25 radii from 0 to 0.25. The
uncorrected centered profile L(d)−d matches the Figure 2 implementation.

Panels H–I compare the frozen Figure 4 realization with Xenium using the
log1p gene-wise raw-count mean and log1p unbiased variance across all 220
genes. Panels J–K report PCC and RMSE between each independent-fit molecular
realization's log1p gene-mean vector and the Xenium vector. Spatial cell-cell
communication inference is not included; co-localization is a noncausal
spatial diagnostic.

## Regeneration

Use the repository environment and SimSpace 0.4.0 (source tag commit
`9889513c0eccd254544a12347c48c0b846e281ba`). From the reproducibility
repository root:

```bash
python supp_figures/SFig9/SFig9.py --prepare-molecular-design
Rscript supp_figures/SFig9/Panel_L_src/generate_molecular_replicates.R \
  main_figures/Fig4/Panel_B_C_D_data/Xenium_reference_metadata.csv \
  main_figures/Fig4/Panel_B_C_D_data/Xenium_reference_count.csv \
  supp_figures/SFig9/Panel_L_data/molecular_simulation_design.tsv \
  supp_figures/SFig9/Panel_L_data/molecular_replicate_summaries.tsv
python supp_figures/SFig9/SFig9.py
```

The R stage was run with R 4.4.2, scDesign3 1.5.0, and
SingleCellExperiment 1.28.1. The final Python stage checks that molecular
summary cell counts and seed provenance match the newly generated spatial
design before plotting.

The workflow writes `supp_figures/SFig9/SFig9.png`, the identical
`example_output/SFig9/SFig9.png`, and frozen result/provenance tables under
`Panel_A_I_data/`, `Panel_J_K_data/`, and `Panel_L_data/`.
`Panel_A_I_data/summary_metrics.tsv` is the numeric source for the manuscript
and response. `analysis_config.json`, `input_manifest.tsv`,
`calibration_manifest.tsv`, `seed_provenance.tsv`, and software-version tables
record the complete regeneration provenance available from the archived fits.
