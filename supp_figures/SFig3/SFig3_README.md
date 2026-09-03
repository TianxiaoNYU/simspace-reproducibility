# Supplementary Figure 3 — finite-basis fixed spatial effects and observation layers

This analysis addresses reviewer comments R2-M1 and R2-m5 with one
coordinated evidence package. It validates SimSpace's optional direct
spatial-expression component and observation controls. It also exports
conditional gene truth that could support a future, separately prespecified
known-truth spatially variable gene (SVG) benchmark. It is not itself a
cross-method SVG benchmark, and the existing R2-M3 concordance/stability
analysis remains separate.

## Scientific question and model

The native reference-free generator first draws phenotype-conditioned gene
means. For selected genes, the optional direct layer now follows the regression
fixed-effect formulation reviewed in *Categorization of 34 computational
methods to detect spatially variable genes from spatially resolved
transcriptomics data* ([Yan, Hua, and Li,
2025](https://doi.org/10.1038/s41467-025-56080-w)). With phenotype 0 configured
as the treatment-coded reference, the full finite-basis
shared-plus-phenotype-deviation validation model is

\[
\log \mu_{ig}=\beta_{0g}+\sum_{k=1}^{3}c_{ik}\beta_{kg}
+B(s_i)^\mathsf{T}\eta_{g0}
+\sum_{k=1}^{3}c_{ik}B(s_i)^\mathsf{T}\eta_{gk}.
\]

Here, \(B(s_i)\) is a simple finite spatial basis, \(\eta_{g0}\) is the
base spatial surface identified by reference phenotype 0, and \(\eta_{gk}\)
is phenotype \(k\)'s deviation from that surface. Thus phenotype 0 has surface
\(\eta_{g0}\), whereas phenotype \(k\) has total surface
\(\eta_{g0}+\eta_{gk}\). The reference interaction is zero for
identifiability. The native library-size factor is one. This is implemented
with the API key `overall_coefficients`, nested gene-to-state
`cell_type_coefficients`, and `reference_state=0`; the compatibility
`coefficients` key is an alias for the base/reference surface in SimSpace.

The validation asks whether the configured base/reference and cell-type-specific
coefficient blocks can be recovered from counts without supplying the generated
phenotype means as an oracle offset, whether likelihood-ratio tests have useful
power and appropriate null rejection rates at unadjusted per-gene
\(\alpha=0.05\), whether observation parameters control
their declared quantities, and whether omitted or zero-valued options preserve
the baseline output exactly.

## Frozen simulation design

- SimSpace distribution: 0.4.2, from source commit
  `747bb234020f807c8fd9963310cd687dd70f1925`. The distribution version, import
  path, and SHA-256 hashes of relevant source files are written to
  `analysis_config.json`.
- Seeds: 0–19; the simulation seed is the independent statistical unit.
- Layout: 50 × 50 reference-free lattice, four phenotypes, one fixed MRF
  affinity matrix, three Gibbs sweeps, 0.65 uniform retention (approximately
  1,600 retained cells per seed), and coordinate jitter standard deviation
  0.55 lattice units.
- Molecular profiles: 500 phenotype-conditioned Gamma–Poisson genes.
- Direct-effect validation: linear, radial, hotspot, and
  distance-to-structure bases are assessed in separate simulations. Each basis
  has 96 nonzero-effect genes and 24 null controls.
- Nonzero allocation: the 96 genes are split evenly among the base/reference
  term and
  deviations for phenotypes 1, 2, and 3. Each of the four roles has three genes
  at each signed strength −1, −0.75, −0.5, −0.25, 0.25, 0.5, 0.75, and 1.
- Null allocation: six genes are assigned to each of the same four test roles.
- Linear coefficients use the fixed unit diagonal direction
  \((1,1)/\sqrt{2}\); vector estimates are projected onto this configured
  direction. The other bases are one-dimensional.
- Capture efficiencies: 1.0, 0.85, 0.7, 0.55, and 0.4.
- Ambient rates: 0, 5, 20, 50, and 100 expected added counts per cell.
- Mean-dependent dropout: off and intercepts −0.5, 0.5, 1.5, and 2.5 with
  slope −1 on `log1p(latent mean)`.
- Combined stress case: capture efficiency 0.7, ambient rate 50, and dropout
  intercept 1.5/slope −1. This condition is not fitted to a real platform.

The baseline Poisson draw is the biological count realization in this model;
capture, ambient background, and dropout are optional subsequent observation
artifacts.

## Fixed-effect estimation and tests

For every validation gene, the analysis fits a Poisson log-link model with an
intercept, indicators for phenotypes 1–3, the complete configured basis
\(B\), and all three phenotype-by-basis interactions. It estimates all
coefficients jointly by maximum likelihood. No generated phenotype-specific
mean is supplied to the fit.

The configured block is the base/reference basis block for an `overall`-API
role gene or the
target interaction block for a cell-type-specific gene. A nested model drops
that complete block while retaining all other terms. Twice the log-likelihood
difference is compared with a chi-square distribution with degrees of freedom
equal to the basis dimension (2 for linear; 1 otherwise). `detected` denotes
an unadjusted per-gene \(p<0.05\). Because the interaction model uses treatment
coding, the base/reference-block LRT tests the spatial surface identified by
phenotype 0 conditional on the non-reference deviation blocks; it is not a
population-average spatial-effect test.

Acceptance checks cover full-rank designs, optimizer convergence, nesting of
the likelihoods, exact recovery from latent means, coefficient bias and sign,
null coefficient bias, the per-gene type-I-error/null rejection rate, power for
\(|\theta|\ge 0.5\), observation-control calibration, and exact no-op
compatibility.

## Archived 20-seed results

All 9,600 full/reduced fixed-effect fits converged and every design was full
rank. The complete latent spatial vector, including off-target zero blocks,
matched the configured vector to a maximum absolute error of
\(1.68\times10^{-14}\). Observed coefficients recovered the configured sign in
99.60% of nonzero gene-level fits; the maximum role-specific median absolute
bias was 0.0595, and the maximum absolute median null estimate was 0.0320.

At the unadjusted per-gene \(\alpha=0.05\), the pooled null rejection rate was
4.95%. Basis-specific rates ranged from 3.13% to 7.29%, and the maximum
basis-by-role rate was 8.33%. Minimum basis-by-role-by-sign power across
\(|\theta|\ge0.5\) was 93.89%. Pooled power was 66.09%, 95.63%, 99.74%, and
100% at absolute coefficients 0.25, 0.5, 0.75, and 1, respectively. These are
component-level, unadjusted test results; they are not FDR-controlled discovery
claims.

Omitted, explicitly disabled, and all-zero direct-effect configurations agreed
exactly across all 20 seeds. Maximum median control errors were 0.00012 for
capture efficiency and 0.0248 ambient counts per cell; maximum seed-level
dropout-probability error was 0.00118.

## Figure panels

- **A:** Seed-0 phenotype labels, the type-1-specific direct log-fold truth,
  and observed Gene 46 counts after combined artifacts. Non-target phenotypes
  are masked in gray in the observed-count map so its color scale represents
  the target phenotype only; gray does not denote a zero count. Gene 46 is
  selected by the display-only, count-blind criterion described below.
- **B:** Separate recovery displays for the base/reference spatial block and pooled
  cell-type-specific deviation blocks across four bases. The lower display
  reports seed-level target-block detection rates by absolute strength; the
  zero level is the null rejection rate and nonzero levels are power.
- **C:** One-factor calibration of capture, ambient background, and dropout.
  Dashed lines are identity lines.
- **D:** Consequences of the strongest one-factor settings and the combined
  condition for library size, zero fraction, the gene-level mean–variance
  association, and latent–observed correlation.

Cells and genes within a simulation are not treated as independent biological
replicates. Coefficients and detection fractions are first summarized within
seed where displayed, and the 20 simulation seeds are the replication units.

For this illustrative panel only, the nine seed-0 linear genes with a configured
cell-type-specific deviation of +1 were ranked without using realized counts.
Six count-blind visualization components were percentile-ranked within this
candidate set and combined: configured truth span (22%), a Poisson working
post-artifact contrast proxy (25%), direct-effect contribution to expected
log-count variation (23%), expected target-cell visibility (12%), a square-root
Poisson working design-information proxy after nuisance projection (10%), and
high-target versus non-target prominence (8%). The two Poisson working proxies
are ranking heuristics rather than exact variance or Fisher-information
calculations for the dropout-contaminated observation law. Gene 46 ranked first
(weighted-percentile score 77.6; configured 90% truth span 0.963; Poisson
working contrast proxy 1.71; direct-effect fraction 0.543; high-target to
non-target 95th-percentile ratio 1.94) and carries a type-1 diagonal interaction.
The weights were fixed for this revision but were not preregistered. This display
selection does not enter the coefficient-recovery or likelihood-ratio summaries.

## Run

From the reproducibility-repository root, use the default environment, which
installs SimSpace 0.4.2:

```bash
conda activate simspace-repro
python supp_figures/SFig3/SFig3.py
```

With no flag, the command reruns the complete 20-seed analysis, evaluates all
acceptance checks, writes the data, and renders the figure. To redraw the
figure from archived tables without rerunning simulations:

```bash
python supp_figures/SFig3/SFig3.py --render-only
```

The installed distribution version is the authoritative release check. The
import path, adjacent source state, and source hashes are also retained for
reproducibility.

## Outputs

- `SFig3.png` — rendered figure;
- `Panel_A_D_data/spatial_effect_recovery.tsv` — one row per seed, basis, and
  validation gene, including configured/latent/estimated vectors, projected
  strengths, full and reduced likelihoods, likelihood-ratio statistics,
  p-values, detections, fit diagnostics, and role labels;
- `Panel_A_D_data/observation_control_calibration.tsv` — raw control/realized
  values for every observation sweep;
- `Panel_A_D_data/observation_metrics.tsv` — seed-level count-property and
  latent-recovery metrics;
- `Panel_A_D_data/backward_compatibility.tsv` — exact default/disabled/zero
  comparisons;
- `Panel_A_D_data/gene_truth.tsv` — base/reference and cell-type-specific spatial
  truth plus phenotype-specific baseline means for every seed and gene;
- `Panel_A_D_data/representative_candidate_screen.tsv` — the complete ranked
  seed-0 display screen for the nine eligible Panel A genes;
- `Panel_A_D_data/representative_map.tsv.gz` — seed-0 values used in Panel A;
- `Panel_A_D_data/acceptance_checks.tsv` — all acceptance criteria;
- `Panel_A_D_data/summary_metrics.tsv` — seed-level coefficient, detection,
  observation, and compatibility summaries; and
- `Panel_A_D_data/analysis_config.json` — complete model, parameters, and
  software/source provenance.

An identical copy of the rendered PNG is written to `example_output/SFig3/`.

## Interpretation limits

This is a controlled component validation, not evidence that one configuration
reproduces the full technical behavior of a particular platform and not a
comparison of SVG methods. The direct layer is opt-in and applies only to the
native reference-free generator; it is not imposed on scDesign3, SRTsim,
Splatter, or other external adapter outputs. The exported truth distinguishes
base/reference, cell-type-specific, null-control, and cell-type-only genes. A broader
cross-method known-truth benchmark would require a separate prespecified
analysis; the existing R2-M3 concordance/stability analysis remains separate.
