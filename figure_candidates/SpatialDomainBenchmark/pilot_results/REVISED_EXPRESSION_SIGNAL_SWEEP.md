# Revised expression-baseline domain-signal sweep

The frozen expression-only PCA(50) plus k-means baseline was evaluated on the revised heterogeneous 50-gene/domain program at seven nominal fold changes for three independent SimSpace layout/count seeds. Within a seed, all random draws are paired across fold changes; only the domain-program multiplier changes.

First tested fold change with ARI >= 0.9: seed 0: 3.5, seed 1: 3.5, seed 2: 4.0.

Revised tier check: hard FC 2.0 median ARI 0.002; moderate FC 3.0 median ARI 0.301 (seed range 0.004--0.381); easy FC 4.0 median ARI 1.000.

## Adjusted Rand index

| Fold change | Seed 0 | Seed 1 | Seed 2 | Median |
|---|---|---|---|---|
| 1.0 | 0.000 | -0.007 | -0.003 | -0.003 |
| 1.5 | 0.002 | 0.002 | -0.004 | 0.002 |
| 2.0 | -0.002 | 0.002 | 0.002 | 0.002 |
| 2.5 | 0.002 | 0.001 | -0.000 | 0.001 |
| 3.0 | 0.004 | 0.381 | 0.301 | 0.301 |
| 3.5 | 1.000 | 1.000 | 0.696 | 1.000 |
| 4.0 | 1.000 | 1.000 | 1.000 | 1.000 |

## Normalized mutual information

| Fold change | Seed 0 | Seed 1 | Seed 2 | Median |
|---|---|---|---|---|
| 1.0 | 0.005 | 0.004 | 0.003 | 0.004 |
| 1.5 | 0.007 | 0.006 | 0.005 | 0.006 |
| 2.0 | 0.005 | 0.005 | 0.008 | 0.005 |
| 2.5 | 0.003 | 0.005 | 0.005 | 0.005 |
| 3.0 | 0.004 | 0.419 | 0.378 | 0.378 |
| 3.5 | 1.000 | 1.000 | 0.720 | 1.000 |
| 4.0 | 1.000 | 1.000 | 1.000 | 1.000 |

This is an engineering sweep over three seeds. The tested grid locates the transition but does not establish a universal detection threshold. The heterogeneous model spreads recovery across seeds and introduces a partial moderate regime, but hard PCA/k-means assignments can still change abruptly within an individual seed.
