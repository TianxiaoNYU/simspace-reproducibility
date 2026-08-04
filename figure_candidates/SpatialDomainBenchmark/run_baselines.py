#!/usr/bin/env python3
"""Run one or both prespecified diagnostic controls on a dataset."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA

from pilot_common import N_DOMAINS, load_dataset, run_with_record, write_assignments


BASELINES = ("expression_pca_kmeans", "coordinates_spectral")


def expression_labels(counts, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Apply the frozen expression-only PCA and k-means baseline."""
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    normalized = counts.multiply(1e4 / np.maximum(library_size, 1.0)[:, None])
    normalized.data = np.log1p(normalized.data)
    embedding = PCA(n_components=50, random_state=seed).fit_transform(
        normalized.toarray()
    )
    labels = KMeans(
        n_clusters=N_DOMAINS, random_state=seed, n_init=20
    ).fit_predict(embedding)
    return labels, embedding


def expression_baseline(
    dataset_dir: Path, output_dir: Path, seed: int
) -> dict[str, object]:
    counts, _, coordinates = load_dataset(dataset_dir)
    labels, embedding = expression_labels(counts, seed)
    write_assignments(
        output_dir, coordinates["location_id"], labels, embedding
    )
    return {
        "software_version": "scikit-learn",
        "clustering": "expression-only PCA(50) + k-means K=4",
        "selected_gene_count": int(counts.shape[1]),
    }


def coordinate_baseline(
    dataset_dir: Path, output_dir: Path, seed: int
) -> dict[str, object]:
    _, _, coordinates = load_dataset(dataset_dir)
    xy = coordinates[["x", "y"]].to_numpy(dtype=float)
    labels = SpectralClustering(
        n_clusters=N_DOMAINS,
        affinity="nearest_neighbors",
        n_neighbors=6,
        assign_labels="kmeans",
        random_state=seed,
        n_init=20,
    ).fit_predict(xy)
    write_assignments(output_dir, coordinates["location_id"], labels, xy)
    return {
        "software_version": "scikit-learn",
        "clustering": "coordinates-only 6-NN spectral clustering, K=4",
        "selected_gene_count": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=100000)
    parser.add_argument(
        "--baselines", nargs="+", choices=BASELINES, default=BASELINES
    )
    args = parser.parse_args()
    callbacks = {
        "expression_pca_kmeans": expression_baseline,
        "coordinates_spectral": coordinate_baseline,
    }
    for method in dict.fromkeys(args.baselines):
        callback = callbacks[method]
        output = args.output_root / method
        run_with_record(
            method,
            args.dataset,
            output,
            args.seed,
            lambda callback=callback, output=output: callback(
                args.dataset, output, args.seed
            ),
        )


if __name__ == "__main__":
    main()
