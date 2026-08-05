# Spatial-domain benchmark figure design

The figure follows the visual grammar of Figure 5 and Supplementary Figure 12: spatial maps, a
composition heatmap, and method-colored line summaries. Panel interiors contain
only subplot titles, axis text, metric labels, and legends. Standalone panels
omit letters; the reproducible complete assembly adds them automatically.

## Panel A — simulated spatial structures

- Curved-layer and irregular-MRF domain maps from the predeclared moderate,
  seed-0 datasets.
- Domain-by-cell-type composition heatmap, averaged across the two displayed
  datasets.
- The two maps occupy an aligned title/map row and the heatmap spans a second
  row, giving Panel A approximately the same height as Panel C.
- Composition colors use a fixed proportion scale from 0 to 0.20.
- The four domain colors are shared with Panel B.

## Panel B — representative domain assignments

- Rows: curved layers and irregular MRF.
- Columns: truth, expression-only control, GraphST, STAGATE, BANKSY, SpaGCN,
  spCLUE, and stLearn.
- The displayed datasets are moderate signal, layout seed 0.
- Predicted cluster IDs are mapped to truth colors by maximum-overlap Hungarian
  assignment for visualization only. Stored predictions and metrics are not
  changed.
- ARI for the displayed dataset is printed directly below every map.

## Panel C — domain and boundary recovery

- A 2 × 2 grid combines the former Panels C and D.
- Columns: curved layers and irregular MRF.
- Rows: adjusted Rand index and one-hop-tolerant boundary F1.
- X-axis: domain-specific expression strength, shown as 2×, 3×, and 4×
  (hard, moderate, easy).
- Lines and markers: median metric value across independent layout seeds.
- Transparent bands: observed minimum-to-maximum seed range.
- One method legend is shared by all four axes.
- Failed or invalid runs are excluded rather than assigned an accuracy of zero.
- Run-validity counts and the incomplete spCLUE run are intentionally omitted
  from the panel and retained only in this note, the workspace README, and the
  aggregate report.

NMI, runtime, memory, and expanded failure details remain in the aggregate
report and can be presented in a supplementary figure or table.

## Complete assembly

- Top row: Panel A at left and Panel C at right.
- Bottom row: the wide Panel B spans the complete figure width.
- The assembled `Fig6.png/.svg` is generated from the same plotting
  primitives as the standalone panels, preserving vector labels in the SVG.
