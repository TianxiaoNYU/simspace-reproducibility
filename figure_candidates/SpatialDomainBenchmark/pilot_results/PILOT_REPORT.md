# Spatial-domain MRF benchmark pilot report

Named-method completion: **6/6 runs**.

This is an engineering and design pilot, not the production benchmark. Metrics are shown to detect broken adapters or trivial simulations; they must not be used to tune method parameters against truth.

## Design diagnostics

- irregular_mrf_seed0_intermediate: intermediate signal with domain fold change 2.0; 68 x 68 grid, 64.9% retention, jitter SD 0.55, and 29.2% truth-boundary nodes in the fixed 6-NN graph; niche MRF phi 5.25, 40 sweeps, Manhattan radius 2; cell-type probability mass 40% enriched / 60% shared, with realized domain/cell-type NMI 0.009.

| Dataset | Method | Status | ARI | NMI | Boundary F1 | Seconds |
|---|---|---|---|---|---|---|
| irregular_mrf_seed0_intermediate | GraphST | success | 0.895 | 0.845 | 0.915 | 88.8 |
| irregular_mrf_seed0_intermediate | STAGATE | success | 0.882 | 0.830 | 0.924 | 229.1 |
| irregular_mrf_seed0_intermediate | BANKSY | success | 0.119 | 0.127 | 0.306 | 16.4 |
| irregular_mrf_seed0_intermediate | SpaGCN | success | 0.974 | 0.955 | 0.998 | 5.4 |
| irregular_mrf_seed0_intermediate | spCLUE | success | 0.855 | 0.799 | 0.909 | 19.7 |
| irregular_mrf_seed0_intermediate | stLearn | success | 0.001 | 0.004 | 0.694 | 5.1 |
| irregular_mrf_seed0_intermediate | Expression baseline | success | 0.002 | 0.007 | 0.700 | 0.6 |
| irregular_mrf_seed0_intermediate | Coordinate baseline | success | 0.154 | 0.155 | 0.270 | 0.3 |

## Failures

No named-method failures.

## Pilot-specific deviations and checks

- BANKSY uses its native mclust endpoint with `lambda=0.8` and `G=4`.
- stLearn uses only physical-distance and expression SME weights; the unused morphology placeholder is checked for output invariance.
- All method metrics are recomputed here from held-out `truth.tsv`; runner input files contain no truth labels.
- The low-signal pilots show that spatially fragmented predictions can receive a high one-hop-tolerant boundary F1 despite near-zero ARI/NMI; boundary F1 must therefore be interpreted jointly with partition metrics.
