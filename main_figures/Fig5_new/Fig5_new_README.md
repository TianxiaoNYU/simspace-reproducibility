# Provisional new Figure 5: spatial-domain benchmarking

This folder contains the compact, reproducible handoff for the spatial-domain
benchmark developed for reviewer comment R2-M7. It is provisionally named
`Fig5_new` and intended to replace the current main-text Figure 5; final
numbering can be resolved later.

The benchmark compares GraphST, STAGATE, BANKSY, SpaGCN, spCLUE, and stLearn
against an expression-only PCA plus k-means negative control. It uses 18
reference-free SimSpace datasets: two spatial patterns, three signal levels,
and three paired layout seeds, with four domains and ten cell types in every
dataset.

## Figure panels

- **Panel A:** spatial structures and domain-specific cell-type compositions.
  The two spatial maps occupy the first row and the composition heatmap spans
  the second row, matching Panel C's height in the assembled figure.
- **Panel B:** truth and method assignments for the predeclared moderate,
  seed-0 representative of each topology; ARI is shown below each map.
- **Panel C:** ARI and boundary-recovery F1 across signal levels and topologies.
  Lines show the median across three layout seeds and bands show the observed
  seed range.

Panel letters are omitted from the standalone panel files. The reproducibly
assembled `Fig5_new` output labels the panels and places A and C across the top
with the wide Panel B spanning the second row. Publication-resolution PNGs and
editable SVGs are retained at the figure-folder root.

## Directory layout

```text
Fig5_new/
  Fig5_new.py                  # single plotting entry point
  Panel_A.png/.svg
  Panel_B.png/.svg
  Panel_C.png/.svg
  Fig5_new.png/.svg            # complete A/C-over-B assembly
  figure_manifest.json         # checksums and plotting provenance
  Panel_A_B_data/
    representative_maps.tsv.gz # compact representative-map cache
  Panel_C_data/
    experiment_summary.tsv     # frozen 126-row metric table
    EXPERIMENT_REPORT.md
    evaluation_manifest.json
    experiment_manifest.json
    experiment_plan.json
  Panel_A_B_C_src/             # full benchmark generation/method pipeline
```

The large count matrices, method embeddings, per-run logs, and other
deterministically regenerated artifacts are excluded. The two compact plotting
inputs above are sufficient to reproduce all three panels in the default
figure environment.

## Render the panels

From the repository root:

```bash
conda activate simspace-repro
python main_figures/Fig5_new/Fig5_new.py
```

The script writes `Panel_A`, `Panel_B`, and `Panel_C` as PNG and SVG files,
assembles `Fig5_new.png/.svg`, and refreshes `figure_manifest.json`. The
complete figure is generated whenever all three panels are requested. A subset
can be rendered with, for example, `--panels B C`.

Predicted cluster identifiers are Hungarian-aligned to truth-domain colors for
Panel B display only. This does not alter stored assignments or metrics.

## Benchmark design

The production matrix is:

- patterns: curved layers and irregular MRF;
- signal levels: hard (2-fold), moderate (3-fold), and easy (4-fold)
  domain-specific expression strength;
- layout seeds: 0, 1, and 2;
- methods: six named spatial methods plus one expression-only negative control.

Within each pattern and layout seed, latent draws are paired across signal
levels. Both patterns use 3,000 locations, 2,000 genes, four domains, ten cell
types, 40% domain-enriched cell-type mass, and 60% shared mass. The molecular
model uses 50 heterogeneous domain-program genes per domain. Full scientific
and implementation details are frozen in
`Panel_A_B_C_src/BENCHMARK_DESIGN.md` and
`Panel_A_B_C_src/benchmark-config.yml`.

The frozen aggregate has all 126 expected rows. One preexisting spCLUE output
for curved-layer, hard, seed 1 contained three rather than four nonempty
clusters and is retained as an invalid row in the aggregate. This validity
detail is documentation only and is not annotated in the main figure.

## Full benchmark reproduction

The production pipeline is kept under `Panel_A_B_C_src/`. It uses two Python
environments, which is the minimum compatible split: one legacy Scanpy stack
for GraphST, STAGATE, SpaGCN, and spCLUE, and one current stack for stLearn.
The generation and evaluation scripts use `simspace-repro`.

Create the shared method environment and run its pinned source/host-R install:

```bash
mamba env create -f main_figures/Fig5_new/Panel_A_B_C_src/environment.yml
conda activate spatial-domain-benchmark
bash main_figures/Fig5_new/Panel_A_B_C_src/install-methods.sh
```

Create the isolated stLearn environment:

```bash
mamba env create -f main_figures/Fig5_new/Panel_A_B_C_src/environment-stlearn.yml
conda activate spatial-domain-benchmark-stlearn
python main_figures/Fig5_new/Panel_A_B_C_src/smoke-test-stlearn.py
```

No Conda R runtime is used. BANKSY and `mclust` use the normal local R library;
the validated host had R 4.4.2, BANKSY 1.2.0, and mclust 6.1.1. Exact
macOS-arm64 environment snapshots are also included as
`resolved-osx-arm64*.yml`.

Inspect the complete command plan without running methods:

```bash
python main_figures/Fig5_new/Panel_A_B_C_src/run_experiment.py --dry-run
```

Run the full production benchmark:

```bash
python main_figures/Fig5_new/Panel_A_B_C_src/run_experiment.py
```

Generated matrices and per-method results are written under
`Panel_A_B_C_src/experiment_data/` and
`Panel_A_B_C_src/experiment_results/`; both paths are gitignored. After a full
rerun, copy the new `experiment_summary.tsv` into `Panel_C_data/`, then rebuild
the representative-map cache and panels with:

```bash
python main_figures/Fig5_new/Fig5_new.py --rebuild-map-cache
```

The no-image stLearn adapter follows the official stSME clustering workflow:
filtering, library normalization, log transformation, PCA, physical-distance
plus expression SME adjustment, scaling, a second PCA, and four-cluster
k-means. BANKSY uses the prespecified `lambda = 0.5` setting.

## Method environments and versions

| Method | Version/source | Immutable revision |
|---|---|---|
| GraphST | 1.1.1 | `d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74` |
| STAGATE_pyG | 1.0.0 | `ae1158ca8cf1eb6bb8ee198298552d44c9ac21db` |
| SpaGCN | 1.2.7 | `dc7a1c26ea0fdf4dfe7064adc7699be141b4871f` |
| spCLUE | source | `bbd2c342e7a67c1617275f721cec2e3f4c23a799` |
| stLearn | 1.4.1 | `03e38844a599ef4f3df8235e3f57bd74a9e96fcf` |
| BANKSY | 1.2.0 | host R library |

`smoke-test.py` validates the four-method environment, compiled PyTorch
Geometric operations, and the local-R bridge. `smoke-test-stlearn.py` validates
the isolated stLearn runtime. Both passed on macOS arm64 during preparation.
