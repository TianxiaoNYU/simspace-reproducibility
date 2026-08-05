# Supplementary Figure 3 — optional spatial-expression and observation layers

This analysis addresses reviewer comments R2-M1 and R2-m5 with one
coordinated evidence package. It also exports conditional gene truth in a
form that can be reused for the related R2-M3 spatially variable gene (SVG)
benchmark. It does not modify or reinterpret the original SimSpace results.

## Scientific question

The original reference-free molecular generator draws phenotype-conditioned
Gamma–Poisson counts. SimSpace 0.3.4 retains that generator unchanged and adds
two independent, opt-in layers:

1. a direct spatial log-mean effect for selected genes,
   \(\mu_{ig}=\mu^{(0)}_{ig}\exp\{B(x_i,y_i)\gamma_g\}\); and
2. an observation layer containing optional binomial capture thinning,
   Poisson ambient background, and excess dropout.

The analysis asks whether the direct coefficient is recoverable from observed
counts, whether each observation parameter controls its declared quantity,
and whether omitted or zero-valued options preserve the historical output
exactly.

## Frozen design

- SimSpace: version 0.3.4, source commit
  `de0a4c002e4ae733e354e3e180ab69b381ad994a`.
- Seeds: 0–19; the simulation seed is the independent statistical unit.
- Layout: 50 × 50 reference-free lattice, four phenotypes, one fixed MRF
  affinity matrix, three Gibbs sweeps, 0.65 uniform retention (approximately
  1,600 retained cells per seed), and coordinate jitter standard deviation
  0.55 lattice units.
- Molecular profiles: 500 phenotype-conditioned Gamma–Poisson genes. The
  linear analysis assigns direct effects to 96 genes and retains 404
  cell-type-only genes. Gene 95 is the representative anti-diagonal-effect
  gene.
- Direct-effect validation: four supported centered, unit-range basis families
  (linear, radial, hotspot, and distance-to-structure). Each basis is assessed
  in its own simulation using 96 nonzero-effect genes and 24 null-control
  genes. Linear effects span horizontal, vertical, diagonal, and anti-diagonal
  directions.
- Direct-effect strengths: −1.0, −0.75, −0.5, −0.25, 0, 0.25, 0.5, 0.75,
  and 1.0. The zero level consists of genes excluded from the opt-in effect.
- Capture efficiencies: 1.0, 0.85, 0.7, 0.55, and 0.4.
- Ambient rates: 0, 5, 20, 50, and 100 expected added counts per cell.
- Mean-dependent dropout: off and intercepts −0.5, 0.5, 1.5, and 2.5 with
  slope −1 on `log1p(latent mean)`.
- Combined stress case: capture efficiency 0.7, ambient rate 50, and dropout
  intercept 1.5/slope −1. This condition was declared before examining the
  results and is not fitted to a real dataset.

The baseline Poisson draw is the biological count realization in this model;
the optional capture, ambient, and dropout components are additional
observation artifacts rather than a relabeling of Poisson variability as
technical noise.

## Figure panels

- **A:** Representative seed-0 phenotype labels, the direct Gene 95 log-fold
  truth, and observed Gene 95 counts after the combined artifact condition.
- **B:** Configured versus recovered direct coefficient across linear, radial,
  hotspot, and structure-distance bases. Each small point is one simulation
  seed after taking the within-seed median across genes assigned to the same
  basis and coefficient; large markers and ranges show the across-seed median
  and full seed range. Open black diamonds show the analytically recovered
  latent coefficient.
- **C:** One-factor calibration of capture, ambient background, and dropout.
  Dashed lines are identity lines.
- **D:** Consequences of the strongest one-factor settings and the combined
  condition for library size, zero fraction, the gene-level mean–variance
  association, and latent–observed correlation.

Cells and genes within one simulation are not treated as biological
replicates. Panel B uses a Poisson log-link fit with the known
phenotype-conditioned baseline mean as an offset; therefore, the recovered
coefficient estimates the direct within-cell-type spatial effect rather than
the marginal association induced by spatially patterned phenotype labels.

Gene 95 was selected without examining its realized spatial count pattern.
Among seed-0 linear genes with configured coefficient +1, it has the largest
minimum phenotype-specific baseline mean; ties are resolved by the median
phenotype-specific baseline mean and then GeneID. This makes expression
visible across all four cell types while avoiding selection on the observed
noise realization.

## Run

From the reproducibility-repository root:

```bash
conda activate simspace-repro
python supp_figures/SFig3/SFig3.py
```

This reruns the complete 20-seed analysis, evaluates all predeclared checks,
and writes both the data and figure. To redraw the figure without rerunning
simulation:

```bash
python supp_figures/SFig3/SFig3.py --render-only
```

## Outputs

- `SFig3.png` — rendered figure; the Python entry script retains the editable
  figure specification;
- `Panel_A_D_data/spatial_effect_recovery.tsv` — one row per seed, basis, and
  validation gene;
- `Panel_A_D_data/observation_control_calibration.tsv` — raw control/realized
  values for every observation sweep;
- `Panel_A_D_data/observation_metrics.tsv` — seed-level count-property and
  latent-recovery metrics;
- `Panel_A_D_data/backward_compatibility.tsv` — exact default/disabled/zero
  comparisons;
- `Panel_A_D_data/gene_truth.tsv` — conditional spatial truth and
  phenotype-specific baseline means for every seed and gene;
- `Panel_A_D_data/representative_map.tsv.gz` — seed-0 values used in Panel A;
- `Panel_A_D_data/acceptance_checks.tsv` — all frozen pass/fail criteria;
- `Panel_A_D_data/summary_metrics.tsv` — canonical seed-level summaries for
  any later manuscript or response numbers; and
- `Panel_A_D_data/analysis_config.json` — complete parameters and software
  provenance.

A copy of the rendered PNG is written to `example_output/SFig3/`.

## Interpretation limits

This is a controlled simulator validation, not a claim that one setting
reproduces the full technical behavior of a particular spatial-transcriptomic
platform. The direct effect is opt-in and applies only to the native
reference-free generator; it is not imposed on scDesign3, SRTsim, Splatter,
or other external adapter outputs. The gene-truth table distinguishes
conditional coordinate effects from cell-type-only genes, but benchmarking
SVG detection methods against that truth remains a separate R2-M3 analysis.
