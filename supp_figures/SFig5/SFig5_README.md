# Supplementary Figure 5 — Gibbs-sweep sensitivity

This analysis addresses reviewer comment R1-8 by varying the
phenotype-level `num_iterations` parameter in the reference-free configuration
used for Supplementary Figure 4. It characterizes the sensitivity of declared
spatial summaries and does not claim a formal convergence diagnostic or proof
of full-field equilibrium.

## Design

- Phenotype Gibbs sweeps: 1, 2, 3, 4, 6, 8, and 10.
- Manuscript default: 4 phenotype sweeps.
- Configurations: eight parameter/seed configurations (seeds 0–7).
- Shared settings: 100 × 100 lattice, two niches, nine cell states,
  Manhattan neighborhood distance 3, six niche-level sweeps, and coordinate
  perturbation standard deviation 0.2.
- Violin-plot summaries: cell-state Moran's I (5-nearest-neighbor graph), Geary's C
  (20-nearest-neighbor graph), and the 50-nearest-neighbor cell-state
  interaction-distance matrix.

For each parameter/seed configuration, the generated parameters and all
settings other than `num_iterations` are held fixed. The script records both the metric
distributions and component-wise agreement with the four-sweep default.
Seeds, cell states, and state pairs are repeated simulation summaries rather
than independent biological replicates.

## Run

The frozen run uses SimSpace 0.4.0. From the reproducibility repository root:

```bash
conda activate simspace-repro
python supp_figures/SFig5/SFig5.py
```

The script verifies the adjacent `../SimSpace` source checkout and writes:

- `supp_figures/SFig5/SFig5.png`;
- `example_output/SFig5/SFig5.png`;
- `Panel_A_B_data/spatial_summary_metrics.tsv`, the frozen raw numeric source;
- distribution and four-sweep agreement summaries; and
- the complete analysis configuration and software versions.

The sensitivity results should be interpreted together with the ten
independently seeded simulations in Supplementary Figure 9, which were
generated from the same archived fitted parameter set.

## Result

Relative to the four-sweep default, the median component-wise Pearson
correlations across the eight configurations at 6, 8, and 10 sweeps were
0.96–0.97 for Moran's I, 0.89–0.96 for Geary's C, and 0.92–0.97 for
interaction distance. The corresponding individual-configuration ranges
were −0.12–0.99, 0.01–0.98, and 0.76–0.98, respectively. Thus, these
declared summaries showed high median agreement around the four-sweep
setting, with configuration-specific variation. This result does not
establish convergence of the complete label field.

The existing figure workflows use four phenotype sweeps and six niche sweeps
for the quantitative 2-D analyses. The representative Supplementary Figure 4C
layouts use six phenotype sweeps (with six niche sweeps), and the Figure 2
3-D example uses five phenotype sweeps and four niche sweeps.
