# Expression-baseline domain-signal sweep

The frozen expression-only PCA(50) plus k-means baseline was evaluated at seven domain fold changes for three independent SimSpace layout/count seeds. Within a seed, all random draws are paired across fold changes; only the domain-program multiplier changes.

First tested fold change with ARI >= 0.9: seed 0: 2.2, seed 1: 2.2, seed 2: >2.2.

## Adjusted Rand index

| Fold change | Seed 0 | Seed 1 | Seed 2 | Median |
|---|---|---|---|---|
| 1.0 | 0.000 | -0.007 | -0.003 | -0.003 |
| 1.2 | 0.004 | -0.006 | 0.001 | 0.001 |
| 1.4 | -0.001 | -0.004 | -0.001 | -0.001 |
| 1.6 | 0.003 | 0.009 | 0.001 | 0.003 |
| 1.8 | 0.001 | 0.000 | 0.002 | 0.001 |
| 2.0 | 0.002 | 0.001 | -0.004 | 0.001 |
| 2.2 | 1.000 | 1.000 | 0.690 | 1.000 |

## Normalized mutual information

| Fold change | Seed 0 | Seed 1 | Seed 2 | Median |
|---|---|---|---|---|
| 1.0 | 0.005 | 0.004 | 0.003 | 0.004 |
| 1.2 | 0.005 | 0.007 | 0.008 | 0.007 |
| 1.4 | 0.005 | 0.004 | 0.005 | 0.005 |
| 1.6 | 0.004 | 0.008 | 0.008 | 0.008 |
| 1.8 | 0.006 | 0.006 | 0.005 | 0.006 |
| 2.0 | 0.007 | 0.005 | 0.005 | 0.005 |
| 2.2 | 1.000 | 1.000 | 0.714 | 1.000 |

This is an engineering sweep over three seeds. The tested grid locates the transition but does not establish a universal detection threshold.
