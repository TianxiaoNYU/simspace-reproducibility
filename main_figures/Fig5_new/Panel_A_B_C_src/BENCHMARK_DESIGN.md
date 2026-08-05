# Spatial-domain benchmark design

Status: the 18-dataset production experiment and all aggregate evaluations are
complete. Standalone main-figure candidate panels are implemented; final figure
numbering remains undecided.

## Scientific question

How accurately and reliably do six spatial-domain methods recover known
spatial domains from fully reference-free SimSpace data containing exactly
2,000 genes?

The estimand is **spatial-domain recovery**. SimSpace niche labels are the
domain truth. Simulated cell phenotypes are a separate biological layer and
must never be used as the evaluation target. This distinction prevents a
cell-type clustering result from being presented as domain detection.

## Primary design

The benchmark contains 18 datasets, 108 named-method runs, and 18
expression-only negative-control runs:

- two pattern types: curved layers and irregular MRF domains;
- three domain-signal levels: easy, moderate, and hard;
- three independent layout seeds per pattern and signal level; and
- six methods per dataset: GraphST, STAGATE, BANKSY, SpaGCN, spCLUE, and
  stLearn.

For each pattern and seed, the coordinates, domain labels, phenotype labels,
baseline molecular means, library-size factors, and observation-noise draws
are paired across the three difficulty levels. Only the multiplier on the
predeclared domain-expression program changes. This makes the difficulty
comparison attributable to domain signal rather than to a new random tissue.

| Property | Frozen primary value |
|---|---:|
| Locations per dataset | 3,000 |
| Genes presented to every method | 2,000 |
| Spatial domains | 4 |
| Cell phenotypes | 10 |
| Domain-enriched phenotypes | 3, 2, 3, and 2 (mean 2.5/domain) |
| Pattern types | curved layers and irregular MRF |
| Layout seeds | 0--2 within each pattern |
| Nominal domain fold change, easy | 4.0 |
| Nominal domain fold change, moderate | 3.0 |
| Nominal domain fold change, hard | 2.0 |
| Named-method runs | 18 datasets x 6 methods = 108 |
| Expression-only control runs | 18 datasets x 1 control = 18 |
| Total evaluated runs | 126 |

The irregular-MRF pattern uses SimSpace's niche-level MRF on a 68 x 68 nominal
grid. Uniform seeded thinning retains exactly 3,000 of 4,624 locations (64.9%),
followed by coordinate jitter with standard deviation 0.55 lattice units,
matching the sampling treatment used for SFig12 panel A. The niche MRF uses
`phi=5.25`, 40 Gibbs sweeps, and a Manhattan radius-2 neighborhood; these
settings were chosen from truth-geometry diagnostics before method rerunning to
replace the initial salt-and-pepper realization with continuous but irregular
domains.

The curved-layer pattern uses four perturbed radial-depth bands in a seeded
reference-free mask on a 104 x 118 nominal grid, retains exactly 3,000
locations, and applies coordinate jitter with standard deviation 0.15 lattice
units. The complete revised molecular model and all three difficulty settings
are identical to those used for the MRF pattern; only the spatial domain
geometry changes. Every domain in either pattern must contain at least 10% of
locations. This is a generation constraint, not post-hoc dataset selection.

## Molecular truth

The 2,000 generic genes are allocated before simulation:

- 1,200 background genes;
- 600 phenotype-program genes, 60 per phenotype; and
- 200 domain-program genes, 50 per domain.

Each cell type has one primary domain under a balanced 3--2--3--2 assignment,
giving ten types and a mean of 2.5 enriched types per domain. The assignment is
rotated by layout seed so that a particular generic type or domain label is not
systematically associated with higher within-domain diversity. Within each
domain, 40% of cell-type probability is divided among its enriched types and
60% is shared across all ten types. Thus, no type is exclusive to a domain and
cell-type labels are not interchangeable with domain labels.

Phenotype-program effects are intrinsic to cell type and remain unchanged when
a type occurs outside its primary domain. The composition structure is held
fixed across the paired easy, moderate, and hard datasets. The revised domain
program scales the log of the nominal fold change by three fixed latent terms:
rotated domain strengths (0.80, 0.95, 1.05, and 1.20), gene-specific weights
from a scaled Beta(2,2) distribution on 0.60--1.40, and location-specific
activity from a scaled Beta(5,2) distribution on 0.50--1.00. These latent
draws are paired across difficulty levels within each layout seed. This reduces
the number of domain genes, breaks the symmetric block effect, and introduces
within-domain heterogeneity while preserving an interpretable nominal fold
change. The same phenotype mixture and domain-expression model is used for
both pattern types. Counts are drawn with SimSpace's reference-free
Gamma--Poisson machinery and one fixed observation-noise configuration. No
real count matrix, real coordinate, real tissue label, reference-derived
marker list, Moran's I target, or other fitted summary enters generation.

The production generator must export the latent mean and gene-class table so
that the molecular signal is auditable. Those files are evaluation provenance
and are not exposed to method runners.

## Information and preprocessing contract

Every method receives the same three items:

1. the raw 3,000 x 2,000 integer count matrix;
2. the two spatial coordinates; and
3. the known domain count, `K = 4`.

No method receives domain truth, phenotype truth, topology name, difficulty
name, image, histology feature, or generative parameter. Giving `K` is an
explicit known-K benchmark choice because several named workflows require a
cluster count. Any resolution search may inspect only the number of produced
clusters, never ARI, NMI, boundaries, or truth labels.

To avoid AnnData-version coupling between the two Python environments, the
canonical interchange format is Matrix Market plus tabular metadata:

- `counts.mtx.gz` -- locations by genes;
- `genes.tsv` -- generic gene identifiers only;
- `coordinates.tsv` -- location identifier, x, and y only;
- `truth.tsv` -- domain and phenotype labels, held outside runner inputs; and
- `manifest.json` -- parameters, seeds, dimensions, and SHA-256 checksums.

The held-out provenance also contains `gene_truth.tsv` and a gzip-compressed
`latent_mean.npy.gz` locations-by-genes float32 array. Neither is exposed to a
method adapter.

Every method starts from all 2,000 genes; the benchmark orchestrator performs
no method-specific pre-screening. A runner may apply documented internal
feature selection, normalization, log transformation, scaling, PCA, or graph
construction as part of the method, but it must record the selected gene count
and identifiers. Coordinates are rescaled once so that the median
nearest-neighbor distance is one. The exact preprocessing operations and
resulting graph size are recorded per run.

## Method adapters

The primary result is each method's documented end-to-end domain workflow,
with only the common known-K and no-histology constraints imposed.

| Method | Runtime | Frozen clustering policy |
|---|---|---|
| GraphST | `spatial-domain-benchmark` | GraphST embedding, then documented mclust with `G=4` |
| STAGATE | `spatial-domain-benchmark` | STAGATE embedding, then documented mclust with `G=4` |
| BANKSY | host R/local library | domain mode with fixed `lambda=0.5`, matching the previous SimSpace figures, then native mclust with `G=4` |
| SpaGCN | `spatial-domain-benchmark` | histology disabled; native resolution search targeting four clusters |
| spCLUE | `spatial-domain-benchmark` | `n_clusters=4`, embedding, then its documented mclust path |
| stLearn | `spatial-domain-benchmark-stlearn` | no image; physical-distance plus expression SME normalization, gene scaling, PCA, then k-means with `K=4` |

The official workflows support the task definitions used here: GraphST is a
spatially informed clustering model; STAGATE learns an expression-and-spatial
embedding for domain identification; SpaGCN directly predicts spatial
domains; spCLUE exposes a domain embedding and clustering; BANKSY's
larger neighborhood weight is its domain-segmentation mode; and stLearn
provides no-image spatial/expression SME features followed by fixed-K
clustering. Exact API calls, source revisions, and every nondefault parameter
belong in one versioned configuration file per adapter.

The current stLearn implementation eagerly reads `X_morphology` while building
all candidate weight matrices, even when the selected
`weights_matrix_pd_gd` path does not use morphology. The no-image adapter will
supply a constant one-column placeholder solely to satisfy that internal code
path, select only physical-distance plus expression weights, and assert in
each run that changing the placeholder leaves the chosen weights byte-identical.
After SME normalization, the adapter follows the official stSME clustering
tutorial by assigning the adjusted matrix to `adata.X`, scaling genes, rerunning
PCA, and clustering the resulting `X_pca` representation.

Expression-only PCA(50) plus k-means with `K=4` is run on every production
dataset as a prespecified negative control. It receives the same count matrix
but does not use spatial coordinates. Its 18 runs are reported separately and
are not included in the 108-run named-method denominator. The coordinates-only
spectral baseline remains an engineering diagnostic and is not part of the
production experiment.

## Reproducibility and resource contract

- Dataset seed and method seed are separate and stored in every result.
- CPU is the default execution target; thread count, hardware model, timeout,
  and memory cap are frozen after the pilot and shared wherever supported.
- Proposed initial limits are one thread, 90 minutes, and 32 GB peak RSS per
  run. A GPU may be adopted only for all neural methods under a new explicit
  resource stratum.
- A failed or timed-out run is retained with its exit status and final log.
  There is no silent replacement with another method or parameter set.
- One retry with the identical seed is allowed only for a documented
  infrastructure interruption, not for a numerical or method failure.
- Each runner writes `assignments.tsv`, `embedding.tsv.gz` when available,
  `run.json`, `stdout.log`, and `stderr.log`.

Generation and aggregation run in `simspace-repro`. The four compatible Python
methods run in `spatial-domain-benchmark`, stLearn runs in
`spatial-domain-benchmark-stlearn`, and BANKSY runs from the existing local R
library. No Conda R runtime is introduced.

## Evaluation

Primary accuracy metrics are adjusted Rand index (ARI) and normalized mutual
information (NMI) against the four held-out SimSpace niche labels. Boundary F1
is a co-primary spatial metric. In a fixed six-nearest-neighbor evaluation
graph, a location is a boundary location when at least one neighbor has a
different label. Precision is the fraction of predicted boundary locations
within one graph hop of truth, recall is the fraction of true boundary
locations within one graph hop of a prediction, and F1 is their harmonic mean.
The graph and tolerance are identical for every method.

Operational outcomes are wall time, peak resident memory, completion status,
and produced cluster count. Accuracy is summarized only for completed runs;
completion rate and a worst-case failure sensitivity are shown alongside it so
that unreliable methods are not favored by complete-case reporting.

The independent unit is the simulated layout seed, not locations or genes.
Plots show all three paired layout observations per pattern and difficulty,
with the median. If uncertainty intervals are displayed, they must be labeled
as descriptive layout-seed bootstrap intervals because only three independent
layouts are available per cell. No single composite accuracy/runtime
leaderboard is computed.

## Completed pilot and freeze gate

Before the 108 production runs, the engineering workflow was required to:

1. generate reference-free curved-layer and irregular-MRF pilot datasets;
2. run all six adapters and both diagnostic baselines;
3. verify exactly 3,000 aligned assignments, four nonempty predicted clusters,
   finite embeddings where applicable, and truth isolation from runner inputs;
4. use pilot timing only to set the common timeout and memory cap;
5. freeze adapter configurations and their checksums; and
6. rerun the pilot from clean environments before launching production.

The pilot cannot be used to choose a method-specific parameter by ARI, NMI, or
visual agreement. If an official default is unusable, the correction must be
justified from the pinned method documentation and applied uniformly to every
dataset before unblinding production results.

### Engineering-pilot outcome (2026-08-03)

The preliminary two-topology engineering rehearsal completed successfully.
Strong expression-only separation in the first curved-layer example was caused
by intentionally coupling layers to distinct cell-type programs. That coupling
has now been replaced in both patterns by the 40% enriched/60% shared
composition and heterogeneous 50-gene/domain program described above. The
curved-layer pattern is therefore restored to the production comparison rather
than discarded on the basis of the obsolete pilot generator.

The initial 55 x 55 MRF grid retained 3,000 of 3,025 positions and used only
0.15 lattice units of coordinate jitter, leaving a visually artificial grid.
The revised pilot increases the grid to 68 x 68, lowers retention to 64.9%, and
raises jitter to 0.55. A second geometry-only calibration increased the niche
MRF segregation from `phi=2.0` and 10 sweeps to `phi=5.25` and 40 sweeps, using
a balanced deterministic niche seed; method results were not used to select
these parameters. Cell-type composition was simultaneously weakened from
70% enriched/30% shared probability mass to 40% enriched/60% shared. The only
adapter repairs retained from the first rehearsal are a GraphST package-layout
import and a `setuptools<81` compatibility pin for the pinned legacy UMAP stack.

The continuous-domain pilot initially used the moderate domain fold change of
2.5, for which the expression-only diagnostic baseline recovered the truth
perfectly. The rerun therefore uses the predeclared hard fold change of 1.6 to
test a nontrivial expression regime without altering any method adapter.
This hard-signal run also exposed that a spatially fragmented partition can
obtain a deceptively high one-hop-tolerant boundary F1 because almost every
location becomes a predicted boundary. Boundary F1 will therefore be reported
only alongside ARI/NMI, and its tolerance definition must be revisited before
freezing the production figure.

Because the paired 2.5 and 1.6 pilots produced sharply different expression
baselines, a third geometry-matched engineering pilot at fold change 2.0 is
used to resolve the transition. Generation and orchestration accept the fold
change and difficulty label explicitly so all three count matrices can be
retained and compared without changing geometry, truth, or method parameters.

A follow-up expression-only sweep tested fold changes 2.0--2.5 in increments
of 0.1 for layout/count seeds 0, 1, and 2, with random draws paired within each
seed. The first tested fold change with expression-baseline ARI at least 0.9
was 2.1, 2.2, and 2.3 for the three seeds, respectively. All seeds had near-zero
ARI at 2.0 and perfect ARI/NMI from 2.3 through 2.5. Thus, the 2.5 result is
reproducible rather than a stale-run artifact, but it lies beyond a sharp,
seed-dependent detection transition. This engineering result should guide the
eventual production signal grid but does not establish a universal threshold.

The sweep was then extended downward to fold changes 1.0--2.2 in increments of
0.2. Across all three seeds, ARI remained within 0.01 of zero from 1.0 through
2.0. At 2.2, seeds 0 and 1 reached ARI 1.0, whereas seed 2 reached ARI 0.69.
This broader grid confirms a flat null regime through 2.0 followed by a narrow,
seed-dependent transition rather than gradual improvement across low signals.

Inspection of the archived PCA embeddings showed that the original uniform
100-gene/domain program created an eigenvalue crossing between phenotype- and
domain-dominated partitions. The primary design is therefore revised to the
heterogeneous 50-gene/domain log-effect model described above. Its three-seed
expression-baseline sweep yielded median ARI 0.002, 0.301, and 1.000 at nominal
fold changes 2.0, 3.0, and 4.0, respectively. These values are retained as the
revised hard, moderate, and easy levels. Moderate recovery was seed-variable
(ARI 0.004--0.381), which is desirable for a replicated benchmark but requires
reporting all layout seeds rather than only a pooled mean. The revised model
spreads the transition across seeds and creates a partial regime; it cannot
remove the discontinuity of hard PCA/k-means assignments within a single seed.

## Candidate figure

A compact main-figure candidate contains:

- **A:** curved-layer and irregular-MRF domain maps plus the corresponding
  domain-by-cell-type composition heatmap;
- **B:** preselected moderate curved-layer and MRF truth maps with the
  expression-only control and six named-method predictions;
- **C:** median ARI and boundary F1 across hard, moderate, and easy signal
  levels with observed seed ranges, arranged as a 2 x 2 metric-by-pattern grid.

The complete figure places A and C across the top and lets the wide Panel B
span the second row. Panel A places the two topology maps above the shared
composition heatmap so its height matches Panel C.

NMI, runtime, memory, all 18 truth maps, and the
complete failure table belong in a supplementary companion if space is tight.
The representative dataset is fixed by seed before method results are viewed.
`FIGURE_PANEL_DESIGN.md` records the final panel-level visual and statistical
contract, and `../Fig5_new.py` writes each standalone panel plus the assembled
complete figure as PNG and SVG.

## Acceptance criteria

- all 18 manifests declare 3,000 locations, 2,000 genes, four valid domains,
  and ten cell phenotypes;
- every method consumes byte-identical counts and coordinates for a dataset;
- no runner input contains either truth label;
- all 108 named-method attempts and all 18 expression-control attempts have
  assignments or an explicit failure record;
- metric recomputation from archived assignments reproduces the summary table;
- figure labels distinguish spatial domains from cell phenotypes; and
- the README, configuration files, environment files, logs, and software
  revisions are sufficient to reproduce every plotted value.
