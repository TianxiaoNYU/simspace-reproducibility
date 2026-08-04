"""Shared I/O, resource, and evaluation utilities for the pilot benchmark."""

from __future__ import annotations

import gzip
import hashlib
import json
import platform
import resource
import time
import traceback
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.neighbors import NearestNeighbors


N_DOMAINS = 4


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_dataset(dataset_dir: Path) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Load only the files available to a method runner."""
    with gzip.open(dataset_dir / "counts.mtx.gz", "rb") as handle:
        counts = mmread(handle).tocsr().astype(np.float32)
    genes = pd.read_csv(dataset_dir / "genes.tsv", sep="\t", dtype=str)
    coordinates = pd.read_csv(dataset_dir / "coordinates.tsv", sep="\t")
    if counts.shape != (len(coordinates), len(genes)):
        raise ValueError(
            f"Count shape {counts.shape} does not match "
            f"{len(coordinates)} locations x {len(genes)} genes."
        )
    if coordinates["location_id"].duplicated().any():
        raise ValueError("Duplicate location identifiers in coordinates.tsv")
    if genes["gene_id"].duplicated().any():
        raise ValueError("Duplicate gene identifiers in genes.tsv")
    return counts, genes, coordinates


def make_adata(dataset_dir: Path):
    import anndata as ad

    counts, genes, coordinates = load_dataset(dataset_dir)
    obs = coordinates.set_index("location_id", drop=False)
    var = genes.set_index("gene_id", drop=False)
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.obsm["spatial"] = coordinates[["x", "y"]].to_numpy(dtype=float)
    return adata


def peak_rss_mb() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if platform.system() == "Darwin":
        return value / (1024.0 * 1024.0)
    return value / 1024.0


def write_assignments(
    output_dir: Path,
    location_ids: np.ndarray | pd.Series | list[str],
    labels: np.ndarray | pd.Series | list[object],
    embedding: np.ndarray | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    location_ids = np.asarray(location_ids).astype(str)
    labels = np.asarray(labels).astype(str)
    if len(location_ids) != len(labels):
        raise ValueError("Assignment length does not match the input locations.")
    pd.DataFrame(
        {"location_id": location_ids, "predicted_domain": labels}
    ).to_csv(output_dir / "assignments.tsv", sep="\t", index=False)
    if embedding is not None:
        embedding = np.asarray(embedding, dtype=float)
        if embedding.shape[0] != len(location_ids):
            raise ValueError("Embedding rows do not match input locations.")
        frame = pd.DataFrame(
            embedding,
            columns=[f"dimension_{index + 1}" for index in range(embedding.shape[1])],
        )
        frame.insert(0, "location_id", location_ids)
        frame.to_csv(
            output_dir / "embedding.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )


def run_with_record(
    method: str,
    dataset_dir: Path,
    output_dir: Path,
    seed: int,
    callback: Callable[[], dict[str, object]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    record: dict[str, object] = {
        "method": method,
        "dataset": dataset_dir.name,
        "method_seed": seed,
        "status": "running",
    }
    try:
        details = callback()
        record.update(details)
        record["status"] = "success"
    except Exception as error:
        record.update(
            {
                "status": "failed",
                "error_type": type(error).__name__,
                "error_message": str(error),
                "traceback": traceback.format_exc(),
            }
        )
        raise
    finally:
        record["wall_time_seconds"] = time.perf_counter() - started
        record["peak_rss_mb"] = peak_rss_mb()
        record["finished_utc"] = pd.Timestamp.now(tz="UTC").isoformat()
        (output_dir / "run.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n"
        )


def boundary_f1(
    coordinates: np.ndarray,
    truth: np.ndarray,
    predicted: np.ndarray,
    n_neighbors: int = 6,
) -> tuple[float, float, float]:
    """One-graph-hop tolerant boundary precision, recall, and F1."""
    graph = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(coordinates)
    neighbors = graph.kneighbors(return_distance=False)[:, 1:]
    true_boundary = boundary_mask(neighbors, truth)
    pred_boundary = boundary_mask(neighbors, predicted)
    true_dilated = true_boundary | np.any(true_boundary[neighbors], axis=1)
    pred_dilated = pred_boundary | np.any(pred_boundary[neighbors], axis=1)
    precision = (
        float(np.mean(true_dilated[pred_boundary])) if pred_boundary.any() else 0.0
    )
    recall = (
        float(np.mean(pred_dilated[true_boundary])) if true_boundary.any() else 0.0
    )
    f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
    return precision, recall, f1


def boundary_mask(neighbors: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Mark nodes whose fixed-neighbor graph crosses a label boundary."""
    labels = np.asarray(labels)
    return np.any(labels[neighbors] != labels[:, None], axis=1)


def boundary_fraction(
    coordinates: np.ndarray, labels: np.ndarray, n_neighbors: int = 6
) -> float:
    """Fraction of nodes adjacent to a different label in a fixed KNN graph."""
    graph = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(coordinates)
    neighbors = graph.kneighbors(return_distance=False)[:, 1:]
    return float(boundary_mask(neighbors, labels).mean())
