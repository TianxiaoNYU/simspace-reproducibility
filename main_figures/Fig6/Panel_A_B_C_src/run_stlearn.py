#!/usr/bin/env python3
"""Run the no-image stLearn SME workflow on one benchmark dataset."""

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
import stlearn as st

from pilot_common import N_DOMAINS, make_adata, run_with_record, write_assignments


def run(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    adata = make_adata(dataset_dir)
    coordinates = adata.obsm["spatial"]
    adata.obs["imagecol"] = coordinates[:, 0]
    adata.obs["imagerow"] = coordinates[:, 1]
    adata.obs["array_col"] = coordinates[:, 0]
    adata.obs["array_row"] = coordinates[:, 1]

    st.pp.filter_genes(adata, min_cells=1)
    st.pp.normalize_total(adata, target_sum=1e4)
    st.pp.log1p(adata)
    components = min(50, adata.n_obs - 1, adata.n_vars - 1)
    st.em.run_pca(
        adata,
        n_comps=components,
        random_state=seed,
    )

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
    sme = np.asarray(adata.obsm["raw_SME_normalized"])
    expression = (
        adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    )
    expression_norm = np.linalg.norm(expression)
    sme_relative_l2_change = float(
        np.linalg.norm(sme - expression) / expression_norm
    )
    physical_neighbors = np.asarray(adata.uns["physical_distance"]).sum(axis=1) - 1

    # Invariance check required by the design: weights_matrix_pd_gd must not
    # depend on the placeholder morphology values. Run this before replacing
    # X_pca with the post-SME PCA, because expression PCA is an intended input
    # to the pd_gd weights.
    selected_weights = adata.uns["weights_matrix_pd_gd"].copy()
    adata_check = adata.copy()
    adata_check.obsm["X_morphology"] = (
        np.arange(adata.n_obs, dtype=float)[:, None] + 1
    )
    from stlearn.spatial.sme._weighting_matrix import weight_matrix

    weight_matrix(adata_check, platform="Visium")
    invariant = np.array_equal(
        selected_weights, adata_check.uns["weights_matrix_pd_gd"]
    )
    if not invariant:
        raise RuntimeError("stLearn no-image placeholder changed pd_gd weights.")

    # Follow the upstream stSME clustering tutorial: make the SME-adjusted
    # matrix the active expression matrix, standardize genes, and then rerun
    # PCA.  Scaling is important here because SME operates on log-normalized
    # expression, whose gene-wise variances otherwise remain highly unequal.
    adata.X = sme.copy()
    st.pp.scale(adata)
    st.em.run_pca(
        adata,
        n_comps=components,
        random_state=seed,
    )
    embedding = np.asarray(adata.obsm["X_pca"])
    st.tl.clustering.kmeans(
        adata,
        n_clusters=N_DOMAINS,
        use_data="X_pca",
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

    return {
        "software_version": "stLearn 1.4.1",
        "clustering": (
            "physical-distance + expression SME, gene scaling, PCA, k-means K=4"
        ),
        "selected_gene_count": int(adata.n_vars),
        "pca_components": int(components),
        "morphology_used": False,
        "placeholder_invariance_check": True,
        "mean_physical_neighbor_count": float(np.mean(physical_neighbors)),
        "sme_relative_l2_change": sme_relative_l2_change,
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
