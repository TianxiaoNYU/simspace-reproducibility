#!/usr/bin/env python3
"""Run GraphST, STAGATE, SpaGCN, or spCLUE on one dataset."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/spatial-domain-pilot-numba")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/spatial-domain-pilot-cache")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import scanpy as sc
import torch
from sklearn.decomposition import PCA

from pilot_common import N_DOMAINS, make_adata, run_with_record, write_assignments


METHODS = ("graphst", "stagate", "spagcn", "spclue")


def graphst(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    # GraphST 1.1.1 ships the model in GraphST/GraphST.py but does not
    # re-export it from the package root.
    from GraphST.GraphST import GraphST as GraphSTModel
    from GraphST.utils import clustering

    adata = make_adata(dataset_dir)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False, max_value=10)
    adata.var["highly_variable"] = True
    model = GraphSTModel(
        adata,
        device=torch.device("cpu"),
        epochs=600,
        random_seed=seed,
        datatype="Stereo",
    )
    result = model.train()
    clustering(
        result,
        n_clusters=N_DOMAINS,
        key="emb",
        method="mclust",
        refinement=False,
    )
    write_assignments(
        output_dir,
        result.obs_names,
        result.obs["domain"].astype(str),
        result.obsm["emb"],
    )
    return {
        "software_version": "GraphST 1.1.1",
        "epochs": 600,
        "clustering": "mclust EEE, G=4",
        "selected_gene_count": int(result.var["highly_variable"].sum()),
        "spatial_graph": "GraphST Stereo KNN, n_neighbors=3",
    }


def stagate(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    import STAGATE_pyG as STAGATE

    adata = make_adata(dataset_dir)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.var["highly_variable"] = True
    STAGATE.Cal_Spatial_Net(adata, model="KNN", k_cutoff=6, verbose=True)
    result = STAGATE.train_STAGATE(
        adata,
        n_epochs=1000,
        random_seed=seed,
        device=torch.device("cpu"),
    )
    result = STAGATE.mclust_R(
        result,
        num_cluster=N_DOMAINS,
        used_obsm="STAGATE",
        random_seed=seed,
    )
    write_assignments(
        output_dir,
        result.obs_names,
        result.obs["mclust"].astype(str),
        result.obsm["STAGATE"],
    )
    return {
        "software_version": "STAGATE_pyG 1.0.0",
        "epochs": 1000,
        "clustering": "mclust EEE, G=4",
        "selected_gene_count": int(result.var["highly_variable"].sum()),
        "spatial_graph": "KNN, k=6",
    }


def spagcn(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    import SpaGCN as spg

    adata = make_adata(dataset_dir)
    spg.prefilter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    coordinates = adata.obsm["spatial"]
    adjacency = spg.calculate_adj_matrix(
        x=coordinates[:, 0], y=coordinates[:, 1], histology=False
    )
    length_scale = spg.search_l(p=0.5, adj=adjacency, start=0.01, end=1000, tol=0.01)
    if length_scale is None:
        raise RuntimeError("SpaGCN could not determine a spatial length scale.")
    np.random.seed(seed)
    torch.manual_seed(seed)
    model = spg.SpaGCN()
    model.set_l(length_scale)
    model.train(
        adata,
        adjacency,
        init="kmeans",
        n_clusters=N_DOMAINS,
        num_pcs=50,
        max_epochs=2000,
    )
    labels, _ = model.predict()
    embedding, _ = model.model.predict(model.embed, model.adj_exp)
    if hasattr(embedding, "detach"):
        embedding = embedding.detach().cpu().numpy()
    write_assignments(output_dir, adata.obs_names, labels, embedding)
    return {
        "software_version": "SpaGCN 1.2.7",
        "max_epochs": 2000,
        "clustering": "native DEC initialized with k-means, K=4",
        "selected_gene_count": int(adata.n_vars),
        "histology": False,
        "length_scale": float(length_scale),
    }


def spclue(dataset_dir: Path, output_dir: Path, seed: int) -> dict[str, object]:
    import spCLUE

    spCLUE.fix_seed(seed)
    adata = make_adata(dataset_dir)
    adata = spCLUE.preprocess(adata)
    components = min(200, adata.n_obs - 1, adata.n_vars - 1)
    adata.obsm["X_pca"] = PCA(
        n_components=components, random_state=seed, svd_solver="randomized"
    ).fit_transform(adata.X)
    graphs = {
        "spatial": spCLUE.prepare_graph(adata, "spatial"),
        "expr": spCLUE.prepare_graph(adata, "expr"),
    }
    model = spCLUE.spCLUE(
        adata.obsm["X_pca"],
        graphs,
        n_clusters=N_DOMAINS,
        epochs=500,
        random_seed=seed,
        device=torch.device("cpu"),
        dim_input=components,
    )
    _, embedding = model.train()
    adata.obsm["spCLUE"] = embedding
    adata = spCLUE.clustering(
        adata,
        N_DOMAINS,
        radius=30,
        key="spCLUE",
        refinement=True,
        cluster_methods="mclust",
    )
    labels = adata.obs["mclust_refined"]
    write_assignments(output_dir, adata.obs_names, labels, embedding)
    return {
        "software_version": "spCLUE commit bbd2c342",
        "epochs": 500,
        "clustering": "mclust EEE, G=4, native 30-neighbor refinement",
        "selected_gene_count": int(adata.n_vars),
        "pca_components": int(components),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=100000)
    args = parser.parse_args()
    callback = globals()[args.method]
    run_with_record(
        args.method,
        args.dataset,
        args.output,
        args.seed,
        lambda: callback(args.dataset, args.output, args.seed),
    )


if __name__ == "__main__":
    main()
