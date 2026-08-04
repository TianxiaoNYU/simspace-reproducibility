"""Import and minimal clustering test for the stLearn benchmark environment."""

from __future__ import annotations

from importlib.metadata import version

import numpy as np
import pandas as pd
import stlearn as st
import torch


rng = np.random.default_rng(7)
spots = [f"spot_{i}" for i in range(8)]
genes = [f"gene_{i}" for i in range(5)]
count = pd.DataFrame(rng.poisson(3, size=(8, 5)), index=spots, columns=genes)
spatial = pd.DataFrame(
    {
        "imagecol": [0, 1, 2, 3, 0, 1, 2, 3],
        "imagerow": [0, 0, 0, 0, 1, 1, 1, 1],
    },
    index=spots,
)

adata = st.create_stlearn(count, spatial, library_id="smoke", scale=1.0)
adata.obsm["X_pca"] = rng.normal(size=(adata.n_obs, 3))
st.tl.clustering.kmeans(
    adata,
    n_clusters=2,
    use_data="X_pca",
    random_state=7,
    key_added="smoke_domain",
)

assert adata.obs["smoke_domain"].nunique() == 2
assert torch.isfinite(torch.ones(2)).all()
print(f"PASS stLearn {version('stlearn')} import and minimal clustering")
print(f"PASS PyTorch {torch.__version__}")
