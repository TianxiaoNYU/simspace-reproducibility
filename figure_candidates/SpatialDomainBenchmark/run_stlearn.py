#!/usr/bin/env python3
"""Run the no-image stLearn SME workflow on one pilot dataset."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/spatial-domain-pilot-numba-stlearn")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/spatial-domain-pilot-cache")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import scanpy as sc
import stlearn as st
from sklearn.decomposition import PCA

from pilot_common import N_DOMAINS, make_adata, run_with_record, write_assignments


def run(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    adata = make_adata(dataset_dir)
    coordinates = adata.obsm["spatial"]
    adata.obs["imagecol"] = coordinates[:, 0]
    adata.obs["imagerow"] = coordinates[:, 1]
    adata.obs["array_col"] = coordinates[:, 0]
    adata.obs["array_row"] = coordinates[:, 1]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    components = min(50, adata.n_obs - 1, adata.n_vars - 1)
    adata.obsm["X_pca"] = PCA(
        n_components=components, random_state=seed
    ).fit_transform(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X)

    # stLearn 1.4.1 builds every candidate weight matrix eagerly, even when
    # weights_matrix_pd_gd is selected. This constant placeholder satisfies
    # the unused morphology branch without conveying image information.
    adata.obsm["X_morphology"] = np.ones((adata.n_obs, 1), dtype=float)
    st.spatial.sme.sme_normalize(
        adata,
        use_data="raw",
        weights="weights_matrix_pd_gd",
        platform="Visium",
    )
    sme = adata.obsm["raw_SME_normalized"]
    embedding = PCA(n_components=components, random_state=seed).fit_transform(sme)
    adata.obsm["X_pca_sme"] = embedding
    st.tl.clustering.kmeans(
        adata,
        n_clusters=N_DOMAINS,
        use_data="X_pca_sme",
        n_init=20,
        random_state=seed,
        key_added="domain",
    )
    write_assignments(
        output_dir,
        adata.obs_names,
        adata.obs["domain"].astype(str),
        embedding,
    )

    # Invariance check required by the design: weights_matrix_pd_gd must not
    # depend on the placeholder morphology values.
    selected_weights = adata.uns["weights_matrix_pd_gd"].copy()
    adata_check = adata.copy()
    adata_check.obsm["X_morphology"] = np.arange(adata.n_obs, dtype=float)[:, None] + 1
    from stlearn.spatial.sme._weighting_matrix import weight_matrix

    weight_matrix(adata_check, platform="Visium")
    invariant = np.array_equal(
        selected_weights, adata_check.uns["weights_matrix_pd_gd"]
    )
    if not invariant:
        raise RuntimeError("stLearn no-image placeholder changed pd_gd weights.")
    return {
        "software_version": "stLearn 1.4.1",
        "clustering": "physical-distance + expression SME, PCA, k-means K=4",
        "selected_gene_count": int(adata.n_vars),
        "pca_components": int(components),
        "morphology_used": False,
        "placeholder_invariance_check": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=100000)
    args = parser.parse_args()
    run_with_record(
        "stlearn",
        args.dataset,
        args.output,
        args.seed,
        lambda: run(args.dataset, args.output, args.seed),
    )


if __name__ == "__main__":
    main()
