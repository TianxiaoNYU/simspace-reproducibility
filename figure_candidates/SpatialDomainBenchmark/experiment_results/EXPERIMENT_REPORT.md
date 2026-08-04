# Spatial-domain production benchmark report

Named-method completion: **107/108 runs** across 18 datasets.

Expression-only negative-control completion: **18/18 runs**.

Each row below summarizes independent layout seeds. Difficulty levels are paired within pattern and seed; only the nominal domain-expression fold change differs within each trio.

## Accuracy and runtime summary

| Pattern | Difficulty | Method | Completed | Median ARI | Median NMI | Median boundary F1 | Median seconds |
|---|---|---|---:|---:|---:|---:|---:|
| curved layers | hard | GraphST | 3/3 | 0.928 | 0.881 | 0.897 | 88.7 |
| curved layers | hard | STAGATE | 3/3 | 0.970 | 0.946 | 0.997 | 225.2 |
| curved layers | hard | BANKSY | 3/3 | 0.313 | 0.389 | 0.363 | 16.2 |
| curved layers | hard | SpaGCN | 3/3 | 0.010 | 0.009 | 0.432 | 5.6 |
| curved layers | hard | spCLUE | 2/3 | 0.949 | 0.918 | 0.987 | 20.1 |
| curved layers | hard | stLearn | 3/3 | -0.000 | 0.013 | 0.433 | 4.7 |
| curved layers | hard | Expression-only control | 3/3 | 0.004 | 0.013 | 0.428 | 0.6 |
| curved layers | moderate | GraphST | 3/3 | 0.973 | 0.950 | 0.998 | 89.4 |
| curved layers | moderate | STAGATE | 3/3 | 0.973 | 0.949 | 0.997 | 226.4 |
| curved layers | moderate | BANKSY | 3/3 | 0.335 | 0.426 | 0.403 | 15.8 |
| curved layers | moderate | SpaGCN | 3/3 | 0.998 | 0.996 | 1.000 | 5.3 |
| curved layers | moderate | spCLUE | 3/3 | 0.954 | 0.923 | 0.990 | 18.6 |
| curved layers | moderate | stLearn | 3/3 | 0.393 | 0.432 | 0.433 | 4.4 |
| curved layers | moderate | Expression-only control | 3/3 | 0.380 | 0.422 | 0.422 | 0.6 |
| curved layers | easy | GraphST | 3/3 | 0.976 | 0.955 | 1.000 | 88.6 |
| curved layers | easy | STAGATE | 3/3 | 0.974 | 0.951 | 0.997 | 225.8 |
| curved layers | easy | BANKSY | 3/3 | 0.310 | 0.430 | 0.413 | 16.6 |
| curved layers | easy | SpaGCN | 3/3 | 0.998 | 0.994 | 1.000 | 5.1 |
| curved layers | easy | spCLUE | 3/3 | 0.954 | 0.923 | 0.990 | 17.6 |
| curved layers | easy | stLearn | 3/3 | 1.000 | 1.000 | 1.000 | 4.3 |
| curved layers | easy | Expression-only control | 3/3 | 1.000 | 1.000 | 1.000 | 0.6 |
| irregular mrf | hard | GraphST | 3/3 | 0.850 | 0.787 | 0.915 | 89.1 |
| irregular mrf | hard | STAGATE | 3/3 | 0.889 | 0.834 | 0.925 | 225.9 |
| irregular mrf | hard | BANKSY | 3/3 | 0.203 | 0.198 | 0.403 | 16.9 |
| irregular mrf | hard | SpaGCN | 3/3 | 0.004 | 0.009 | 0.679 | 5.2 |
| irregular mrf | hard | spCLUE | 3/3 | 0.850 | 0.793 | 0.903 | 27.2 |
| irregular mrf | hard | stLearn | 3/3 | -0.005 | 0.005 | 0.678 | 4.2 |
| irregular mrf | hard | Expression-only control | 3/3 | 0.002 | 0.007 | 0.681 | 0.6 |
| irregular mrf | moderate | GraphST | 3/3 | 0.905 | 0.851 | 0.926 | 90.1 |
| irregular mrf | moderate | STAGATE | 3/3 | 0.890 | 0.836 | 0.924 | 226.1 |
| irregular mrf | moderate | BANKSY | 3/3 | 0.181 | 0.213 | 0.381 | 17.6 |
| irregular mrf | moderate | SpaGCN | 3/3 | 0.977 | 0.956 | 0.979 | 5.0 |
| irregular mrf | moderate | spCLUE | 3/3 | 0.854 | 0.798 | 0.909 | 17.6 |
| irregular mrf | moderate | stLearn | 3/3 | 0.289 | 0.369 | 0.693 | 4.2 |
| irregular mrf | moderate | Expression-only control | 3/3 | 0.301 | 0.380 | 0.693 | 0.6 |
| irregular mrf | easy | GraphST | 3/3 | 0.906 | 0.856 | 0.921 | 89.2 |
| irregular mrf | easy | STAGATE | 3/3 | 0.891 | 0.835 | 0.922 | 224.7 |
| irregular mrf | easy | BANKSY | 3/3 | 0.179 | 0.236 | 0.392 | 17.1 |
| irregular mrf | easy | SpaGCN | 3/3 | 0.982 | 0.965 | 0.991 | 4.8 |
| irregular mrf | easy | spCLUE | 3/3 | 0.854 | 0.798 | 0.912 | 17.7 |
| irregular mrf | easy | stLearn | 3/3 | 1.000 | 1.000 | 1.000 | 4.2 |
| irregular mrf | easy | Expression-only control | 3/3 | 1.000 | 1.000 | 1.000 | 0.6 |

## Failures and invalid outputs

| Dataset | Method | Status | Error |
|---|---|---|---|
| curved_layers_seed1_hard | spCLUE | invalid_output | ValueError: Expected four nonempty clusters, found 3 |

## Evaluation contract

- ARI and NMI compare predictions with held-out SimSpace domain labels.
- Boundary F1 uses the fixed six-nearest-neighbor, one-hop-tolerant definition.
- A completed result must contain exactly 3,000 aligned assignments and four nonempty clusters.
- Accuracy metrics are omitted for failures and invalid outputs; completion remains explicit.
- `truth.tsv` is read only during this aggregation step, never by a method runner.
- The expression-only control receives counts but no spatial coordinates in its computation and is reported separately from the six named methods.
