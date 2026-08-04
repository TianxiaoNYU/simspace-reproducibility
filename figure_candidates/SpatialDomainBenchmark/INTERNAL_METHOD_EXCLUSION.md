# Internal method-exclusion note

This file records a rebuttal decision and is not manuscript, figure, or caption
text.

## DeepGFT

DeepGFT was attempted during the first engineering pilot but is excluded from
the revised benchmark before the MRF rerun. The decision is based on execution
and reproducibility concerns, not on selecting methods by pilot accuracy:

- the upstream fast Fourier path requests a full eigendecomposition at this
  dataset size; making the pilot runnable required a nondefault bounded sparse
  eigensolver with 256 modes;
- the two 3,000-location pilot runs required 394.4 and 378.7 seconds on CPU,
  making it the slowest adapter in the rehearsal; and
- stable reproducibility has not been reliable enough to justify presenting
  the output as a comparable final benchmark result.

Policy: do not run this method in the revised pilot or production workflow, and
do not include it in manuscript method lists, figures, captions, aggregate
tables, or the public benchmark README. The reviewer response should state
frankly that it was evaluated at the engineering stage but omitted because a
stable, reproducible, runtime-feasible result could not be obtained.

Suggested response language:

> We also evaluated DeepGFT during development. However, we could not obtain a
> sufficiently stable reproducible result, and its Fourier/training
> workflow was substantially more computationally demanding at the benchmark
> scale. We therefore did not include DeepGFT in the final comparison rather
> than report a potentially misleading result from an unstable execution.
