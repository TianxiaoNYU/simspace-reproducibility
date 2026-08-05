# Proposed new Figure 2 (`Fig2_new`)

This directory contains the standalone, code-generated panels for a proposed
main figure between the current Figures 2 and 3. The repository's existing
figure numbering and `FIGURE_INDEX.md` are intentionally unchanged.

## Figure brief

**Takeaway:** reference-free SimSpace parameters can produce a designed,
layered cortex; a controllable family of related spatial tissues; and
cell-type-conditional molecular profiles with biologically motivated cortical
marker patterns.

The spatial generator does not ingest observed coordinates, spatial labels, a
tissue mask, or fitted spatial statistics. The molecular generator does not
ingest an expression count matrix. The named marker profiles are qualitative
biological design parameters, not estimates fitted from mouse scRNA-seq.

## Panels

### Panel A: niche design

`Panel_A.png` shows the four curved niches in the selected anchor realization:

1. superficial (L2/3);
2. middle (L4/5);
3. deep (L5);
4. inner (L5/6-L6b).

The anchor uses the previously selected `104 x 118` lattice and contains 3,360
retained cells. Its complete spatial parameter set is written to
`output/cortex_parameters.json`.

### Panel B: two-parameter cortex family

`Panel_B.png` is a 5 x 5 sweep. All stochastic seeds and all unshown
parameters are held fixed so that the displayed differences are attributable
to the two axes.

- Columns change the `inner_radius` parameter from `0.58` to `0.30`, producing
  progressively thicker cortical sectors.
- Rows change a layer-specificity factor from `3.0` to `0.0`. A factor of `1.0`
  is the selected niche-by-cell-type mixture; `0.0` gives every niche the same
  global cell-type profile; values above `1.0` sharpen niche enrichment.
- The center tile (`inner_radius = 0.44`, specificity `1.0`) is the exact anchor
  used in Panels A and C.

`Panel_B_data/sweep_parameters.tsv` records the plotted parameter combination,
cell count, and informal tuning diagnostics for every tile.

### Panel C: reference-free molecular profile

`Panel_C.png` displays simulated counts for two contrasting cortical markers
on the fixed anchor:

- *Cux2*: L2/3 IT-enriched superficial expression;
- *Foxp2*: L6 CT-enriched deep expression.

The expression design contains 99 genes: nine named cortical subclass markers,
two additional generic markers per cell type, and 72 background genes. SimSpace
draws Poisson counts from the declared cell-type-conditional means; the mean
parameters are generated from Gamma distributions. Spatial expression arises
through the simulated cell labels, without adding a direct coordinate effect.

The nine named marker-to-subclass assignments are:

| Marker | Designed enrichment |
|---|---|
| *Cux2* | L2/3 IT |
| *Rorb* | L4/5 IT |
| *Fezf2* | L5 ET |
| *Deptor* | L5 IT |
| *Tshz2* | L5/6 NP |
| *Foxp2* | L6 CT |
| *Sulf2* | L6 IT |
| *Car3* | L6 IT Car3 |
| *Ctgf* | L6b |

These choices are biologically motivated by mouse cortical cell-taxonomy work,
including [Tasic et al., Nature 2018](https://doi.org/10.1038/s41586-018-0654-5),
[Zhang et al., Nature 2021](https://doi.org/10.1038/s41586-021-03223-w), and the
[BICCN multimodal census, Nature 2021](https://doi.org/10.1038/s41586-021-03950-0).
They are used as illustrative enrichment patterns rather than quantitative
reference-derived values.

## Reproduce the panels

From the `simspace-reproducibility` repository root:

```bash
conda activate simspace-repro
python main_figures/Fig2_new/Fig2_new.py
```

The script deterministically regenerates the three panel PNGs and all compact
panel data. Use `--dpi` to change the Panel B and C raster resolution; Panel A
uses the established 300-dpi renderer from `cortex_reference_free.py`.

## Contents

- `Fig2_new.py`: entry script for Panels A-C
- `cortex_reference_free.py`: de-novo cortex geometry and SimSpace generator
- `Panel_A.png`: four-niche anchor design
- `Panel_B.png`: 5 x 5 thickness-by-specificity sweep
- `Panel_C.png`: *Cux2* and *Foxp2* spatial expression
- `Panel_A_C_data/anchor_cells.tsv.gz`: anchor coordinates, labels, niches, and depth
- `Panel_B_data/sweep_parameters.tsv`: complete 25-tile sweep table
- `Panel_C_data/gene_parameters.tsv`: cell-type-conditional expression means
- `Panel_C_data/simulated_expression.tsv.gz`: 3,360 x 99 simulated count matrix
- `Panel_C_data/marker_expression_by_cell_type.tsv`: marker-expression summary
- `figure_brief.json`: frozen figure claims, axes, assumptions, and anchor parameters
- `output/`: archived tuning panels, anchor files, and visual-comparison diagnostics
