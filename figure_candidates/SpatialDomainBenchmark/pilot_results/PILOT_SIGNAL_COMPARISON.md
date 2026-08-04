# Paired pilot signal comparison

All three pilots use byte-identical coordinates, truth labels, gene identifiers, and gene annotations. Only the domain-expression fold change and resulting count matrix differ. This is a single-layout engineering comparison, not the production benchmark.

## Engineering observations

- The expression baseline remains null at FC 1.6 (ARI 0.003) and FC 2.0 (0.002), then reaches ARI 1.000 at FC 2.5.
- SpaGCN changes earlier, from ARI 0.002 at FC 1.6 to 0.974 at FC 2.0.
- stLearn remains null at FC 2.0 (ARI 0.001) and reaches 1.000 at FC 2.5.
- GraphST, STAGATE, and spCLUE remain strong and comparatively stable across all three paired count matrices; BANKSY remains weak.
- These apparent thresholds come from one layout seed and require replication before biological or method-level interpretation.

## Adjusted Rand index

| Method | FC 1.6 | FC 2.0 | FC 2.5 |
|---|---|---|---|
| GraphST | 0.885 | 0.895 | 0.897 |
| STAGATE | 0.872 | 0.882 | 0.881 |
| BANKSY | 0.116 | 0.119 | 0.177 |
| SpaGCN | 0.002 | 0.974 | 0.979 |
| spCLUE | 0.852 | 0.855 | 0.858 |
| stLearn | 0.005 | 0.001 | 1.000 |
| Expression baseline | 0.003 | 0.002 | 1.000 |
| Coordinate baseline | 0.154 | 0.154 | 0.154 |

## Normalized mutual information

| Method | FC 1.6 | FC 2.0 | FC 2.5 |
|---|---|---|---|
| GraphST | 0.834 | 0.845 | 0.848 |
| STAGATE | 0.819 | 0.830 | 0.829 |
| BANKSY | 0.112 | 0.127 | 0.178 |
| SpaGCN | 0.004 | 0.955 | 0.963 |
| spCLUE | 0.796 | 0.799 | 0.802 |
| stLearn | 0.004 | 0.004 | 1.000 |
| Expression baseline | 0.004 | 0.007 | 1.000 |
| Coordinate baseline | 0.155 | 0.155 | 0.155 |
