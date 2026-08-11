"""Supplementary Figure 9: local reference-guided validation for R1-2/R1-3.

The analysis deliberately reuses the frozen Figure 4 Xenium tile and fitted
SimSpace parameters.  It does not refit the model or treat cells, genes, or
simulation seeds as independent biological replicates.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import sys
from importlib.metadata import version as distribution_version
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/simspace-r1-2-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/simspace-r1-2-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import pearsonr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import simspace as ss


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
FIG3_DATA = REPO_ROOT / "main_figures" / "Fig4" / "Panel_B_C_D_data"
DATA_DIR = SCRIPT_DIR / "Panel_A_I_data"
EXPRESSION_DATA_DIR = SCRIPT_DIR / "Panel_J_K_data"
MOLECULAR_DATA_DIR = SCRIPT_DIR / "Panel_L_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig9"

K_NEIGHBORS = 20
WHOLE_LAYOUT_RIPLEY_RADII = np.linspace(0.0, 0.25, 25)
SIMULATION_SEEDS = tuple(range(10))
PERMUTATIONS = 100
MIN_CELLTYPE_COUNT = 20
MIN_GENE_DETECTION = 0.05
MAX_GENE_DETECTION = 0.95
SIMSPACE_EXPECTED_VERSION = "0.3.4"
SIMSPACE_SOURCE_COMMIT = "de0a4c002e4ae733e354e3e180ab69b381ad994a"

TYPE_ABBREVIATIONS = {
    "Invasive_Tumor": "Tumor",
    "Stromal": "Stromal",
    "Prolif_Invasive_Tumor": "Prolif.",
    "Endothelial": "Endoth.",
    "Macrophages_1": "Macro.",
    "CD8+_T_Cells": "CD8 T",
    "T_Cell_&_Tumor_Hybrid": "T–tumor",
    "B_Cells": "B",
    "CD4+_T_Cells": "CD4 T",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def coordinate_hash(coordinates: np.ndarray) -> str:
    rounded = np.round(np.asarray(coordinates, dtype="<f8"), 8)
    return hashlib.sha256(rounded.tobytes()).hexdigest()


def normalize_coordinates(coordinates: pd.DataFrame | np.ndarray) -> np.ndarray:
    values = np.asarray(coordinates, dtype=float)
    minimum = values.min(axis=0)
    span = np.ptp(values, axis=0)
    span[span == 0] = 1.0
    return (values - minimum) / span


def knn_indices(coordinates: np.ndarray, k: int = K_NEIGHBORS) -> np.ndarray:
    n_cells = len(coordinates)
    if n_cells <= 1:
        raise ValueError("At least two coordinates are required.")
    effective_k = min(k, n_cells - 1)
    indices = cKDTree(coordinates).query(
        coordinates, k=effective_k + 1
    )[1]
    if indices.ndim == 1:
        indices = indices[:, None]
    return indices[:, 1:]


def graph_edges(neighbors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source = np.repeat(np.arange(neighbors.shape[0]), neighbors.shape[1])
    target = neighbors.reshape(-1)
    return source, target


def moran_geary(
    values: np.ndarray, edges: tuple[np.ndarray, np.ndarray]
) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    denominator = float(centered @ centered)
    if denominator <= 0:
        return np.nan, np.nan
    source, target = edges
    n_cells = len(values)
    total_weight = len(source)
    moran_i = (
        n_cells
        / total_weight
        * np.sum(centered[source] * centered[target])
        / denominator
    )
    geary_c = (
        (n_cells - 1)
        / (2 * total_weight)
        * np.sum((values[source] - values[target]) ** 2)
        / denominator
    )
    return float(moran_i), float(geary_c)


def whole_layout_centered_ripley_l(
    coordinates: np.ndarray,
    radii: np.ndarray = WHOLE_LAYOUT_RIPLEY_RADII,
) -> np.ndarray:
    """Figure 2 uncorrected L(r)-r for all points in a normalized unit square."""
    coordinates = np.asarray(coordinates, dtype=float)
    n_points = len(coordinates)
    if n_points < 2:
        return np.full(len(radii), np.nan)
    pairs = cKDTree(coordinates).query_pairs(
        float(radii.max()), output_type="ndarray"
    )
    if len(pairs) == 0:
        return -radii.copy()
    pair_distances = np.linalg.norm(
        coordinates[pairs[:, 0]] - coordinates[pairs[:, 1]], axis=1
    )
    pair_distances.sort()
    pair_counts = np.searchsorted(pair_distances, radii, side="right")
    k_values = 2.0 * pair_counts / (n_points * (n_points - 1))
    return np.sqrt(k_values / np.pi) - radii


def celltype_metrics(
    coordinates: np.ndarray,
    labels: np.ndarray,
    cell_types: list[str],
    dataset: str,
    seed: int | str,
) -> pd.DataFrame:
    neighbors = knn_indices(coordinates)
    edges = graph_edges(neighbors)
    rows: list[dict[str, object]] = []
    for cell_type in cell_types:
        indicator = (labels == cell_type).astype(float)
        moran_i, geary_c = moran_geary(indicator, edges)
        rows.append(
            {
                "dataset": dataset,
                "seed": seed,
                "cell_type": cell_type,
                "n_cells": int(indicator.sum()),
                "moran_i": moran_i,
                "geary_c": geary_c,
            }
        )
    return pd.DataFrame(rows)


def metric_agreement(
    reference: pd.DataFrame,
    simulated: pd.DataFrame,
    seed: int,
    id_column: str,
    metrics: list[str],
) -> pd.DataFrame:
    merged = reference.merge(
        simulated,
        on=id_column,
        suffixes=("_reference", "_simulated"),
        validate="one_to_one",
    )
    rows = []
    for metric in metrics:
        reference_values = merged[f"{metric}_reference"].to_numpy(dtype=float)
        simulated_values = merged[f"{metric}_simulated"].to_numpy(dtype=float)
        valid = np.isfinite(reference_values) & np.isfinite(simulated_values)
        rows.append(
            {
                "seed": seed,
                "metric": metric,
                "pearson_r": float(
                    pearsonr(
                        reference_values[valid], simulated_values[valid]
                    ).statistic
                ),
                "rmse": float(
                    np.sqrt(
                        np.mean(
                            (
                                reference_values[valid]
                                - simulated_values[valid]
                            )
                            ** 2
                        )
                    )
                ),
                "n_units": int(valid.sum()),
            }
        )
    return pd.DataFrame(rows)


def ripley_profile_table(
    dataset: str, seed: int | str, profile: np.ndarray
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "dataset": dataset,
            "seed": seed,
            "normalized_radius": WHOLE_LAYOUT_RIPLEY_RADII,
            "centered_l": profile,
        }
    )


def ripley_profile_agreement(
    reference_profile: np.ndarray, simulated_profile: np.ndarray, seed: int
) -> dict[str, float | int]:
    return {
        "seed": seed,
        "pearson_r": float(
            pearsonr(reference_profile, simulated_profile).statistic
        ),
        "rmse": float(
            np.sqrt(np.mean((reference_profile - simulated_profile) ** 2))
        ),
        "n_radii": len(WHOLE_LAYOUT_RIPLEY_RADII),
    }


def gene_metrics(
    reference_coordinates: np.ndarray,
    simulated_coordinates: np.ndarray,
    reference_counts: pd.DataFrame,
    simulated_counts: pd.DataFrame,
) -> pd.DataFrame:
    reference_log = np.log1p(reference_counts.astype(float))
    simulated_log = np.log1p(simulated_counts.astype(float))
    reference_detection = (reference_counts > 0).mean(axis=0)
    simulated_detection = (simulated_counts > 0).mean(axis=0)
    genes = reference_detection[
        (reference_detection >= MIN_GENE_DETECTION)
        & (reference_detection <= MAX_GENE_DETECTION)
    ].index.intersection(simulated_counts.columns)

    reference_edges = graph_edges(knn_indices(reference_coordinates))
    simulated_edges = graph_edges(knn_indices(simulated_coordinates))
    rows = []
    for gene in genes:
        ref_values = reference_log[gene].to_numpy()
        sim_values = simulated_log[gene].to_numpy()
        ref_moran, ref_geary = moran_geary(ref_values, reference_edges)
        sim_moran, sim_geary = moran_geary(sim_values, simulated_edges)
        rows.append(
            {
                "gene": gene,
                "reference_detection_fraction": float(
                    reference_detection[gene]
                ),
                "simulated_detection_fraction": float(
                    simulated_detection[gene]
                ),
                "moran_i_reference": ref_moran,
                "moran_i_simulated": sim_moran,
                "geary_c_reference": ref_geary,
                "geary_c_simulated": sim_geary,
            }
        )
    return pd.DataFrame(rows)


def gene_agreement(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ["moran_i", "geary_c"]:
        reference = metrics[f"{metric}_reference"].to_numpy()
        simulated = metrics[f"{metric}_simulated"].to_numpy()
        valid = np.isfinite(reference) & np.isfinite(simulated)
        rows.append(
            {
                "metric": metric,
                "pearson_r": float(
                    pearsonr(reference[valid], simulated[valid]).statistic
                ),
                "rmse": float(
                    np.sqrt(
                        np.mean((reference[valid] - simulated[valid]) ** 2)
                    )
                ),
                "n_genes": int(valid.sum()),
            }
        )
    return pd.DataFrame(rows)


def expression_fidelity(
    reference_counts: pd.DataFrame,
    simulated_counts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare gene-wise count means and variances on a log1p scale."""
    genes = reference_counts.columns.intersection(
        simulated_counts.columns, sort=False
    )
    reference = reference_counts.loc[:, genes].astype(float)
    simulated = simulated_counts.loc[:, genes].astype(float)
    metrics = pd.DataFrame(
        {
            "gene": genes,
            "reference_mean_count": reference.mean(axis=0).to_numpy(),
            "simulated_mean_count": simulated.mean(axis=0).to_numpy(),
            "reference_variance": reference.var(axis=0, ddof=1).to_numpy(),
            "simulated_variance": simulated.var(axis=0, ddof=1).to_numpy(),
            "reference_zero_fraction": (reference == 0).mean(axis=0).to_numpy(),
            "simulated_zero_fraction": (simulated == 0).mean(axis=0).to_numpy(),
        }
    )
    for statistic in ["mean_count", "variance"]:
        metrics[f"reference_log1p_{statistic}"] = np.log1p(
            metrics[f"reference_{statistic}"]
        )
        metrics[f"simulated_log1p_{statistic}"] = np.log1p(
            metrics[f"simulated_{statistic}"]
        )

    rows = []
    for statistic in ["mean_count", "variance"]:
        reference_values = metrics[f"reference_log1p_{statistic}"].to_numpy()
        simulated_values = metrics[f"simulated_log1p_{statistic}"].to_numpy()
        valid = np.isfinite(reference_values) & np.isfinite(simulated_values)
        rows.append(
            {
                "statistic": statistic,
                "scale": f"log1p gene-wise {statistic.replace('_', ' ')}",
                "pearson_r": float(
                    pearsonr(
                        reference_values[valid], simulated_values[valid]
                    ).statistic
                ),
                "rmse": float(
                    np.sqrt(
                        np.mean(
                            (
                                reference_values[valid]
                                - simulated_values[valid]
                            )
                            ** 2
                        )
                    )
                ),
                "n_genes": int(valid.sum()),
            }
        )
    return metrics, pd.DataFrame(rows)


def molecular_replicate_agreement(
    reference_counts: pd.DataFrame,
    replicate_summaries: pd.DataFrame,
) -> pd.DataFrame:
    """Compare each molecular realization's gene means with the reference."""
    reference_means = np.log1p(reference_counts.astype(float).mean(axis=0))
    rows = []
    for seed, frame in replicate_summaries.groupby("molecular_seed"):
        simulated_means = frame.set_index("gene")["mean_count"]
        genes = reference_means.index.intersection(
            simulated_means.index, sort=False
        )
        reference_values = reference_means.loc[genes].to_numpy()
        simulated_values = np.log1p(
            simulated_means.loc[genes].to_numpy(dtype=float)
        )
        valid = np.isfinite(reference_values) & np.isfinite(simulated_values)
        rows.append(
            {
                "molecular_seed": int(seed),
                "pearson_r": float(
                    pearsonr(
                        reference_values[valid], simulated_values[valid]
                    ).statistic
                ),
                "rmse": float(
                    np.sqrt(
                        np.mean(
                            (
                                reference_values[valid]
                                - simulated_values[valid]
                            )
                            ** 2
                        )
                    )
                ),
                "n_genes": int(valid.sum()),
            }
        )
    agreement = pd.DataFrame(rows).sort_values("molecular_seed")
    if agreement["molecular_seed"].tolist() != list(SIMULATION_SEEDS):
        raise ValueError(
            "Molecular summaries must contain exactly seeds 0--9."
        )
    return agreement


def centered_colocalization(
    labels: np.ndarray,
    neighbors: np.ndarray,
    cell_types: list[str],
    rng: np.random.Generator,
) -> np.ndarray:
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    one_hot = np.zeros((len(labels), len(cell_types)), dtype=float)
    for row, label in enumerate(labels):
        if label in type_index:
            one_hot[row, type_index[label]] = 1.0

    local_composition = one_hot[neighbors].mean(axis=1)
    observed = local_composition.T @ local_composition / len(labels)
    permutation_mean = np.zeros_like(observed)
    for _ in range(PERMUTATIONS):
        permuted = one_hot[rng.permutation(len(labels))]
        permuted_composition = permuted[neighbors].mean(axis=1)
        permutation_mean += (
            permuted_composition.T
            @ permuted_composition
            / len(labels)
            / PERMUTATIONS
        )
    return observed - permutation_mean


def matrix_long(
    matrix: np.ndarray,
    cell_types: list[str],
    dataset: str,
    seed: int | str,
) -> pd.DataFrame:
    rows = []
    for i, first in enumerate(cell_types):
        for j, second in enumerate(cell_types):
            rows.append(
                {
                    "dataset": dataset,
                    "seed": seed,
                    "cell_type_a": first,
                    "cell_type_b": second,
                    "centered_colocalization": float(matrix[i, j]),
                }
            )
    return pd.DataFrame(rows)


def colocalization_agreement(
    reference: np.ndarray, simulated: np.ndarray, seed: int
) -> dict[str, float | int]:
    upper = np.triu_indices_from(reference, k=1)
    reference_values = reference[upper]
    simulated_values = simulated[upper]
    return {
        "seed": seed,
        "pearson_r": float(
            pearsonr(reference_values, simulated_values).statistic
        ),
        "rmse": float(
            np.sqrt(np.mean((reference_values - simulated_values) ** 2))
        ),
        "n_celltype_pairs": len(reference_values),
    }


def niche_metrics(
    reference_domains: pd.DataFrame,
    simulated_domains: pd.DataFrame,
    cell_types: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ari = adjusted_rand_score(
        simulated_domains["niche"], simulated_domains["Tumor_label"]
    )
    nmi = normalized_mutual_info_score(
        simulated_domains["niche"], simulated_domains["Tumor_label"]
    )
    metric_rows = [
        {
            "metric": "encoded_niche_vs_BANKSY_ARI",
            "domain": "all_simulated_cells",
            "value": float(ari),
            "n_units": len(simulated_domains),
        },
        {
            "metric": "encoded_niche_vs_BANKSY_NMI",
            "domain": "all_simulated_cells",
            "value": float(nmi),
            "n_units": len(simulated_domains),
        },
    ]
    composition_rows = []
    for domain in ["Tumor", "Stromal"]:
        reference = (
            reference_domains.loc[
                reference_domains["Tumor_label"] == domain, "fitted_celltype"
            ]
            .value_counts(normalize=True)
            .reindex(cell_types, fill_value=0.0)
        )
        simulated = (
            simulated_domains.loc[
                simulated_domains["Tumor_label"] == domain,
                "fitted_celltype",
            ]
            .value_counts(normalize=True)
            .reindex(cell_types, fill_value=0.0)
        )
        correlation = float(pearsonr(reference, simulated).statistic)
        rmse = float(np.sqrt(np.mean((reference - simulated) ** 2)))
        metric_rows.extend(
            [
                {
                    "metric": "domain_composition_pearson_r",
                    "domain": domain,
                    "value": correlation,
                    "n_units": len(cell_types),
                },
                {
                    "metric": "domain_composition_rmse",
                    "domain": domain,
                    "value": rmse,
                    "n_units": len(cell_types),
                },
            ]
        )
        for cell_type in cell_types:
            composition_rows.append(
                {
                    "domain": domain,
                    "cell_type": cell_type,
                    "reference_fraction": float(reference[cell_type]),
                    "simulated_fraction": float(simulated[cell_type]),
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(composition_rows)


def summary_rows(
    reference_metadata: pd.DataFrame,
    reference_counts: pd.DataFrame,
    seed_provenance: pd.DataFrame,
    celltype_agreement_table: pd.DataFrame,
    ripley_agreement_table: pd.DataFrame,
    gene_agreement_table: pd.DataFrame,
    expression_agreement_table: pd.DataFrame,
    molecular_replicate_agreement_table: pd.DataFrame,
    colocalization_agreement_table: pd.DataFrame,
    niche_metric_table: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add(metric: str, value: object, unit: str, source: str) -> None:
        rows.append(
            {
                "metric": metric,
                "value": value,
                "unit": unit,
                "source_table": source,
            }
        )

    add("reference_cells", len(reference_metadata), "cells", "input")
    add("reference_genes", reference_counts.shape[1], "genes", "input")
    add(
        "postfit_spatial_realizations",
        len(seed_provenance),
        "seeds",
        "seed_provenance.tsv",
    )
    add(
        "all_coordinate_hashes_unique",
        int(seed_provenance["coordinate_hash"].nunique() == len(seed_provenance)),
        "boolean",
        "seed_provenance.tsv",
    )
    add(
        "maximum_coordinate_matches_to_reference",
        int(seed_provenance["coordinate_matches_to_reference"].max()),
        "coordinates",
        "seed_provenance.tsv",
    )
    for metric in ["moran_i", "geary_c"]:
        values = celltype_agreement_table.loc[
            celltype_agreement_table["metric"] == metric, "pearson_r"
        ]
        add(
            f"celltype_{metric}_pearson_r_median",
            float(values.median()),
            "Pearson r",
            "celltype_agreement.tsv",
        )
        add(
            f"celltype_{metric}_pearson_r_range",
            f"{values.min():.6f},{values.max():.6f}",
            "min,max",
            "celltype_agreement.tsv",
        )
    for column in ["pearson_r", "rmse"]:
        values = ripley_agreement_table[column]
        add(
            f"whole_layout_ripley_{column}_median",
            float(values.median()),
            column,
            "whole_layout_ripley_agreement.tsv",
        )
        add(
            f"whole_layout_ripley_{column}_range",
            f"{values.min():.6f},{values.max():.6f}",
            "min,max",
            "whole_layout_ripley_agreement.tsv",
        )
    for _, row in gene_agreement_table.iterrows():
        add(
            f"gene_{row['metric']}_pearson_r",
            float(row["pearson_r"]),
            "Pearson r",
            "gene_agreement.tsv",
        )
        add(
            f"gene_{row['metric']}_rmse",
            float(row["rmse"]),
            "metric units",
            "gene_agreement.tsv",
        )
    add(
        "eligible_genes",
        int(gene_agreement_table["n_genes"].min()),
        "genes",
        "gene_agreement.tsv",
    )
    for _, row in expression_agreement_table.iterrows():
        add(
            f"expression_log1p_{row['statistic']}_pearson_r",
            float(row["pearson_r"]),
            "Pearson r",
            "../Panel_J_K_data/expression_agreement.tsv",
        )
        add(
            f"expression_log1p_{row['statistic']}_rmse",
            float(row["rmse"]),
            "log1p summary units",
            "../Panel_J_K_data/expression_agreement.tsv",
        )
    add(
        "expression_shared_genes",
        int(expression_agreement_table["n_genes"].min()),
        "genes",
        "../Panel_J_K_data/expression_agreement.tsv",
    )
    for column, unit in [("pearson_r", "Pearson r"), ("rmse", "log1p mean-count units")]:
        values = molecular_replicate_agreement_table[column]
        add(
            f"molecular_replicate_{column}_median",
            float(values.median()),
            unit,
            "../Panel_L_data/molecular_replicate_agreement.tsv",
        )
        add(
            f"molecular_replicate_{column}_range",
            f"{values.min():.6f},{values.max():.6f}",
            "min,max",
            "../Panel_L_data/molecular_replicate_agreement.tsv",
        )
    add(
        "molecular_replicates",
        int(len(molecular_replicate_agreement_table)),
        "seeds",
        "../Panel_L_data/molecular_replicate_agreement.tsv",
    )
    for column in ["pearson_r", "rmse"]:
        values = colocalization_agreement_table[column]
        add(
            f"colocalization_{column}_median",
            float(values.median()),
            column,
            "colocalization_agreement.tsv",
        )
        add(
            f"colocalization_{column}_range",
            f"{values.min():.6f},{values.max():.6f}",
            "min,max",
            "colocalization_agreement.tsv",
        )
    for _, row in niche_metric_table.iterrows():
        add(
            f"niche_{row['metric']}_{row['domain']}",
            float(row["value"]),
            row["metric"],
            "niche_concordance.tsv",
        )
    return pd.DataFrame(rows)


def add_panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )


def plot_figure(
    cell_types: list[str],
    celltype_agreement_table: pd.DataFrame,
    ripley_profile_table_data: pd.DataFrame,
    ripley_agreement_table: pd.DataFrame,
    gene_metric_table: pd.DataFrame,
    gene_agreement_table: pd.DataFrame,
    expression_metric_table: pd.DataFrame,
    expression_agreement_table: pd.DataFrame,
    molecular_replicate_agreement_table: pd.DataFrame,
    reference_colocalization: np.ndarray,
    simulated_colocalization_median: np.ndarray,
    niche_composition_table: pd.DataFrame,
) -> plt.Figure:
    sns.set_theme(style="whitegrid", context="paper")
    fig = plt.figure(figsize=(16.8, 12.6), constrained_layout=True)
    grid = fig.add_gridspec(3, 4)
    bottom_grid = grid[2, :].subgridspec(1, 3)
    spatial_agreement_axis = fig.add_subplot(grid[0, 0])
    gene_axes = [
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
    ]
    ripley_axis = fig.add_subplot(grid[0, 3])
    colocalization_axes = [
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    niche_axis = fig.add_subplot(grid[1, 2])
    expression_axes = [
        fig.add_subplot(grid[1, 3]),
        fig.add_subplot(bottom_grid[0, 0]),
    ]
    robustness_axes = [
        fig.add_subplot(bottom_grid[0, 1]),
        fig.add_subplot(bottom_grid[0, 2]),
    ]

    agreement_plot = pd.concat(
        [
            celltype_agreement_table,
            ripley_agreement_table.assign(metric="whole_layout_ripley"),
        ],
        ignore_index=True,
    )
    agreement_plot["metric_label"] = agreement_plot["metric"].map(
        {
            "moran_i": "Moran's I\ncell types\n(calibration)",
            "geary_c": "Geary's C\ncell types\n(diagnostic)",
            "whole_layout_ripley": "Ripley's L–d\nall cells\n(diagnostic)",
        }
    )
    order = [
        "Moran's I\ncell types\n(calibration)",
        "Geary's C\ncell types\n(diagnostic)",
        "Ripley's L–d\nall cells\n(diagnostic)",
    ]
    sns.boxplot(
        data=agreement_plot,
        x="metric_label",
        y="pearson_r",
        order=order,
        color="#d8e6f3",
        fliersize=0,
        width=0.55,
        ax=spatial_agreement_axis,
    )
    sns.stripplot(
        data=agreement_plot,
        x="metric_label",
        y="pearson_r",
        order=order,
        color="#2f6f9f",
        size=4,
        jitter=0.16,
        ax=spatial_agreement_axis,
    )
    spatial_agreement_axis.axhline(0, color="0.45", linewidth=0.8)
    spatial_agreement_axis.set_ylim(-1.05, 1.05)
    spatial_agreement_axis.set_xlabel("")
    spatial_agreement_axis.set_ylabel("Reference–simulation Pearson r")
    spatial_agreement_axis.set_title("SimSpace vs Xenium (10 seeds)")
    add_panel_label(spatial_agreement_axis, "A")

    gene_labels = {
        "moran_i": ("Moran's I", "B"),
        "geary_c": ("Geary's C", "C"),
    }
    for axis, metric in zip(gene_axes, gene_labels):
        reference = gene_metric_table[f"{metric}_reference"]
        simulated = gene_metric_table[f"{metric}_simulated"]
        axis.scatter(
            reference,
            simulated,
            s=12,
            alpha=0.55,
            color="#2f75b5",
            edgecolors="none",
        )
        low = min(reference.min(), simulated.min())
        high = max(reference.max(), simulated.max())
        padding = 0.04 * (high - low if high > low else 1.0)
        axis.plot(
            [low - padding, high + padding],
            [low - padding, high + padding],
            linestyle="--",
            color="0.4",
            linewidth=0.8,
        )
        agreement = gene_agreement_table.loc[
            gene_agreement_table["metric"] == metric
        ].iloc[0]
        axis.text(
            0.04,
            0.95,
            f"r = {agreement['pearson_r']:.2f}\nRMSE = {agreement['rmse']:.3f}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
        label, panel = gene_labels[metric]
        axis.set_title(f"Gene {label} (n={len(gene_metric_table)})")
        axis.set_xlabel("Xenium reference")
        axis.set_ylabel("SimSpace seed 0")
        add_panel_label(axis, panel)

    reference_ripley = ripley_profile_table_data.loc[
        ripley_profile_table_data["dataset"] == "Xenium_reference"
    ].sort_values("normalized_radius")
    simulated_ripley = ripley_profile_table_data.loc[
        ripley_profile_table_data["dataset"] == "SimSpace"
    ].pivot(
        index="normalized_radius",
        columns="seed",
        values="centered_l",
    )
    radii = simulated_ripley.index.to_numpy(dtype=float)
    ripley_axis.fill_between(
        radii,
        simulated_ripley.min(axis=1).to_numpy(),
        simulated_ripley.max(axis=1).to_numpy(),
        color="#8fb9dc",
        alpha=0.3,
        label="SimSpace range",
    )
    ripley_axis.plot(
        radii,
        simulated_ripley.median(axis=1).to_numpy(),
        color="#2f75b5",
        linewidth=1.8,
        label="SimSpace median",
    )
    ripley_axis.plot(
        reference_ripley["normalized_radius"],
        reference_ripley["centered_l"],
        color="0.2",
        linewidth=1.8,
        label="Xenium reference",
    )
    ripley_axis.axhline(0, color="0.55", linewidth=0.7)
    ripley_axis.text(
        0.04,
        0.95,
        (
            f"median r = {ripley_agreement_table['pearson_r'].median():.2f}\n"
            f"median RMSE = {ripley_agreement_table['rmse'].median():.3f}"
        ),
        transform=ripley_axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
    )
    ripley_axis.set_title("Whole-layout Ripley's L profile")
    ripley_axis.set_xlabel("Normalized radius d")
    ripley_axis.set_ylabel("L(d) − d")
    ripley_axis.legend(frameon=False, fontsize=7, loc="lower left")
    add_panel_label(ripley_axis, "D")

    abbreviations = [TYPE_ABBREVIATIONS.get(t, t) for t in cell_types]
    common_limit = float(
        max(
            np.abs(reference_colocalization).max(),
            np.abs(simulated_colocalization_median).max(),
        )
    )
    for axis, matrix, title, panel in [
        (
            colocalization_axes[0],
            reference_colocalization,
            "Xenium centered co-localization",
            "E",
        ),
        (
            colocalization_axes[1],
            simulated_colocalization_median,
            "SimSpace median centered co-localization",
            "F",
        ),
    ]:
        sns.heatmap(
            matrix,
            cmap="vlag",
            center=0,
            vmin=-common_limit,
            vmax=common_limit,
            square=True,
            xticklabels=abbreviations,
            yticklabels=abbreviations,
            cbar=panel == "F",
            cbar_kws={"label": "Observed − permutation mean", "shrink": 0.7},
            ax=axis,
        )
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=55, labelsize=7)
        axis.tick_params(axis="y", rotation=0, labelsize=7)
        add_panel_label(axis, panel)

    domain_palette = {"Tumor": "#b34f4f", "Stromal": "#3f7f9f"}
    for domain, frame in niche_composition_table.groupby("domain"):
        niche_axis.scatter(
            frame["reference_fraction"],
            frame["simulated_fraction"],
            s=28,
            alpha=0.8,
            color=domain_palette[domain],
            label=domain,
        )
        correlation = pearsonr(
            frame["reference_fraction"], frame["simulated_fraction"]
        ).statistic
        niche_axis.text(
            0.04,
            0.94 if domain == "Tumor" else 0.84,
            f"{domain}: r = {correlation:.3f}",
            color=domain_palette[domain],
            transform=niche_axis.transAxes,
            va="top",
            fontsize=8,
        )
    maximum = max(
        niche_composition_table["reference_fraction"].max(),
        niche_composition_table["simulated_fraction"].max(),
    )
    niche_axis.plot(
        [0, maximum], [0, maximum], linestyle="--", color="0.4", linewidth=0.8
    )
    niche_axis.set_xlabel("Xenium cell-type fraction")
    niche_axis.set_ylabel("SimSpace cell-type fraction")
    niche_axis.set_title("BANKSY domain composition")
    niche_axis.legend(frameon=False, loc="lower right", fontsize=8)
    add_panel_label(niche_axis, "G")

    expression_panels = [
        (
            expression_axes[0],
            "mean_count",
            "Gene-wise mean expression",
            "H",
        ),
        (
            expression_axes[1],
            "variance",
            "Gene-wise expression variance",
            "I",
        ),
    ]
    for axis, statistic, title, panel in expression_panels:
        reference = expression_metric_table[
            f"reference_log1p_{statistic}"
        ]
        simulated = expression_metric_table[
            f"simulated_log1p_{statistic}"
        ]
        axis.scatter(
            reference,
            simulated,
            s=15,
            alpha=0.58,
            color="#2f75b5",
            edgecolors="none",
        )
        low = min(reference.min(), simulated.min())
        high = max(reference.max(), simulated.max())
        padding = 0.04 * (high - low if high > low else 1.0)
        axis.plot(
            [low - padding, high + padding],
            [low - padding, high + padding],
            linestyle="--",
            color="0.4",
            linewidth=0.8,
        )
        agreement = expression_agreement_table.loc[
            expression_agreement_table["statistic"] == statistic
        ].iloc[0]
        axis.text(
            0.04,
            0.95,
            (
                f"PCC = {agreement['pearson_r']:.3f}\n"
                f"RMSE = {agreement['rmse']:.3f}"
            ),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        axis.set_title(
            f"{title} (n={int(agreement['n_genes'])})"
        )
        axis.set_xlabel("Xenium reference (log1p)")
        axis.set_ylabel("SimSpace seed 0 (log1p)")
        add_panel_label(axis, panel)

    seeds = molecular_replicate_agreement_table["molecular_seed"].to_numpy()
    robustness_specs = [
        (
            robustness_axes[0],
            "pearson_r",
            "PCC across 10 molecular seeds",
            "PCC",
            "#2f75b5",
            "J",
            (0.96, 0.05, "right", "bottom"),
        ),
        (
            robustness_axes[1],
            "rmse",
            "RMSE across 10 molecular seeds",
            "RMSE (log1p mean-count units)",
            "#c45a3c",
            "K",
            (0.96, 0.95, "right", "top"),
        ),
    ]
    for axis, metric, title, ylabel, color, panel, annotation in robustness_specs:
        values = molecular_replicate_agreement_table[metric]
        axis.scatter(
            seeds,
            values,
            s=34,
            color=color,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )
        axis.set_xticks(seeds)
        axis.set_xlabel("Simulation seed")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(False)
        text_x, text_y, horizontal_alignment, vertical_alignment = annotation
        axis.text(
            text_x,
            text_y,
            (
                f"median = {values.median():.3f}\n"
                f"range = {values.min():.3f}–{values.max():.3f}"
            ),
            transform=axis.transAxes,
            ha=horizontal_alignment,
            va=vertical_alignment,
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8},
        )
        add_panel_label(axis, panel)

    return fig


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXPRESSION_DATA_DIR.mkdir(parents=True, exist_ok=True)
    MOLECULAR_DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)

    simspace_version = distribution_version("simspace")
    if simspace_version != SIMSPACE_EXPECTED_VERSION:
        raise RuntimeError(
            "Supplementary Figure 9 requires SimSpace "
            f"{SIMSPACE_EXPECTED_VERSION}; found {simspace_version}."
        )

    input_paths = {
        "reference_metadata": FIG3_DATA / "Xenium_reference_metadata.csv",
        "reference_counts": FIG3_DATA / "Xenium_reference_count.csv",
        "simulated_metadata": FIG3_DATA / "simspace_fitted_metadata.csv",
        "simulated_counts": FIG3_DATA / "simspace_fitted_count.csv",
        "fitted_parameters": FIG3_DATA / "simspace_fitted_params.json",
        "reference_banksy_domains": FIG3_DATA / "BANKSY_xenium_domain.csv",
        "simulated_banksy_domains": FIG3_DATA / "BANKSY_simspace_domain.csv",
    }
    manifest = pd.DataFrame(
        [
            {
                "role": role,
                "path_from_repository_root": str(path.relative_to(REPO_ROOT)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for role, path in input_paths.items()
        ]
    )
    manifest.to_csv(DATA_DIR / "input_manifest.tsv", sep="\t", index=False)

    reference_metadata = pd.read_csv(input_paths["reference_metadata"])
    reference_counts = pd.read_csv(
        input_paths["reference_counts"], index_col=0
    ).T
    reference_counts = reference_counts.loc[
        reference_metadata["cell_id"].astype(str)
    ]
    frozen_simulated_metadata = pd.read_csv(
        input_paths["simulated_metadata"], index_col=0
    )
    simulated_counts = pd.read_csv(
        input_paths["simulated_counts"], index_col=0
    )
    if len(frozen_simulated_metadata) != len(simulated_counts):
        raise ValueError("Figure 4 simulated metadata/count rows do not align.")

    state_to_celltype = (
        frozen_simulated_metadata[["state", "fitted_celltype"]]
        .drop_duplicates()
        .set_index("state")["fitted_celltype"]
        .to_dict()
    )
    if len(state_to_celltype) != frozen_simulated_metadata["state"].nunique():
        raise ValueError("A simulated state maps to more than one cell type.")

    reference_coordinates = normalize_coordinates(
        reference_metadata[["x_centroid", "y_centroid"]]
    )
    reference_labels = reference_metadata["Cluster"].to_numpy()
    reference_coordinate_set = {
        tuple(row)
        for row in np.round(reference_coordinates * 100.0, 6)
    }
    cell_types = (
        reference_metadata["Cluster"]
        .value_counts()
        .loc[lambda counts: counts >= MIN_CELLTYPE_COUNT]
        .index.tolist()
    )

    reference_celltype = celltype_metrics(
        reference_coordinates,
        reference_labels,
        cell_types,
        dataset="Xenium_reference",
        seed="reference",
    )
    reference_ripley_profile = whole_layout_centered_ripley_l(
        reference_coordinates
    )
    reference_neighbors = knn_indices(reference_coordinates)
    reference_colocalization = centered_colocalization(
        reference_labels,
        reference_neighbors,
        cell_types,
        np.random.default_rng(9000),
    )

    celltype_tables = [reference_celltype]
    agreement_tables = []
    ripley_profile_tables = [
        ripley_profile_table(
            "Xenium_reference", "reference", reference_ripley_profile
        )
    ]
    ripley_agreements = []
    colocalization_tables = [
        matrix_long(
            reference_colocalization,
            cell_types,
            dataset="Xenium_reference",
            seed="reference",
        )
    ]
    colocalization_agreements = []
    simulated_matrices = []
    seed_rows = []
    molecular_design_tables = []
    seed_zero_metadata: pd.DataFrame | None = None

    for seed in SIMULATION_SEEDS:
        simulation = ss.util.sim_from_json(
            input_file=str(input_paths["fitted_parameters"]),
            shape=(100, 100),
            num_iteration=4,
            n_iter=6,
            seed=seed,
        )
        metadata = simulation.meta.copy()
        metadata.columns = ["state", "row", "col", "niche", "state_rank"]
        metadata["state"] = metadata["state"].astype(int)
        metadata["fitted_celltype"] = metadata["state"].map(state_to_celltype)
        if metadata["fitted_celltype"].isna().any():
            raise ValueError(f"Unmapped cell state in simulation seed {seed}.")
        coordinates = normalize_coordinates(metadata[["row", "col"]])
        labels = metadata["fitted_celltype"].to_numpy()
        molecular_design = metadata[
            ["row", "col", "fitted_celltype"]
        ].copy()
        molecular_design.insert(0, "cell_index", np.arange(len(metadata)))
        molecular_design.insert(0, "molecular_seed", seed)
        molecular_design_tables.append(molecular_design)

        if seed == 0:
            same_state = np.array_equal(
                metadata["state"].to_numpy(),
                frozen_simulated_metadata["state"].astype(int).to_numpy(),
            )
            same_coordinates = np.allclose(
                metadata[["row", "col"]].to_numpy(),
                frozen_simulated_metadata[["row", "col"]].to_numpy(),
            )
            if not (same_state and same_coordinates):
                raise ValueError(
                    "Seed-0 replay does not match the frozen Figure 4 metadata."
                )
            seed_zero_metadata = metadata.copy()

        coordinate_matches = len(
            reference_coordinate_set.intersection(
                {
                    tuple(row)
                    for row in np.round(coordinates * 100.0, 6)
                }
            )
        )
        seed_rows.append(
            {
                "simulation_seed": seed,
                "n_cells": len(metadata),
                "coordinate_hash": coordinate_hash(coordinates),
                "coordinate_matches_to_reference": coordinate_matches,
                "parameter_file_sha256": manifest.loc[
                    manifest["role"] == "fitted_parameters", "sha256"
                ].iloc[0],
            }
        )

        simulated_ripley_profile = whole_layout_centered_ripley_l(coordinates)
        ripley_profile_tables.append(
            ripley_profile_table("SimSpace", seed, simulated_ripley_profile)
        )
        ripley_agreements.append(
            ripley_profile_agreement(
                reference_ripley_profile, simulated_ripley_profile, seed
            )
        )

        metrics = celltype_metrics(
            coordinates,
            labels,
            cell_types,
            dataset="SimSpace",
            seed=seed,
        )
        celltype_tables.append(metrics)
        agreement_tables.append(
            metric_agreement(
                reference_celltype[["cell_type", "moran_i", "geary_c"]],
                metrics[["cell_type", "moran_i", "geary_c"]],
                seed=seed,
                id_column="cell_type",
                metrics=["moran_i", "geary_c"],
            )
        )

        simulated_colocalization = centered_colocalization(
            labels,
            knn_indices(coordinates),
            cell_types,
            np.random.default_rng(10_000 + seed),
        )
        simulated_matrices.append(simulated_colocalization)
        colocalization_tables.append(
            matrix_long(
                simulated_colocalization,
                cell_types,
                dataset="SimSpace",
                seed=seed,
            )
        )
        colocalization_agreements.append(
            colocalization_agreement(
                reference_colocalization, simulated_colocalization, seed
            )
        )

    if seed_zero_metadata is None:
        raise RuntimeError("Seed 0 was not generated.")

    seed_provenance = pd.DataFrame(seed_rows)
    celltype_metric_table = pd.concat(celltype_tables, ignore_index=True)
    celltype_agreement_table = pd.concat(
        agreement_tables, ignore_index=True
    )
    whole_layout_ripley_table = pd.concat(
        ripley_profile_tables, ignore_index=True
    )
    whole_layout_ripley_agreement_table = pd.DataFrame(ripley_agreements)
    colocalization_table = pd.concat(
        colocalization_tables, ignore_index=True
    )
    colocalization_agreement_table = pd.DataFrame(
        colocalization_agreements
    )
    molecular_design_table = pd.concat(
        molecular_design_tables, ignore_index=True
    )
    molecular_design_table.to_csv(
        MOLECULAR_DATA_DIR / "molecular_simulation_design.tsv",
        sep="\t",
        index=False,
    )

    simulated_coordinates = normalize_coordinates(
        seed_zero_metadata[["row", "col"]]
    )
    simulated_labels = seed_zero_metadata["fitted_celltype"].to_numpy()
    gene_metric_table = gene_metrics(
        reference_coordinates,
        simulated_coordinates,
        reference_counts,
        simulated_counts,
    )
    gene_agreement_table = gene_agreement(gene_metric_table)
    expression_metric_table, expression_agreement_table = expression_fidelity(
        reference_counts,
        simulated_counts,
    )
    molecular_summary_path = (
        MOLECULAR_DATA_DIR / "molecular_replicate_summaries.tsv"
    )
    if not molecular_summary_path.exists():
        raise FileNotFoundError(
            f"Missing {molecular_summary_path}. Run "
            "Panel_L_src/generate_molecular_replicates.R first."
        )
    molecular_replicate_summaries = pd.read_csv(
        molecular_summary_path, sep="\t"
    )
    molecular_replicate_agreement_table = molecular_replicate_agreement(
        reference_counts,
        molecular_replicate_summaries,
    )

    reference_domains = pd.read_csv(input_paths["reference_banksy_domains"])
    simulated_domains = pd.read_csv(input_paths["simulated_banksy_domains"])
    niche_metric_table, niche_composition_table = niche_metrics(
        reference_domains,
        simulated_domains,
        sorted(reference_metadata["Cluster"].unique()),
    )

    summary = summary_rows(
        reference_metadata,
        reference_counts,
        seed_provenance,
        celltype_agreement_table,
        whole_layout_ripley_agreement_table,
        gene_agreement_table,
        expression_agreement_table,
        molecular_replicate_agreement_table,
        colocalization_agreement_table,
        niche_metric_table,
    )

    seed_provenance.to_csv(
        DATA_DIR / "seed_provenance.tsv", sep="\t", index=False
    )
    celltype_metric_table.to_csv(
        DATA_DIR / "celltype_spatial_metrics.tsv", sep="\t", index=False
    )
    celltype_agreement_table.to_csv(
        DATA_DIR / "celltype_agreement.tsv", sep="\t", index=False
    )
    whole_layout_ripley_table.to_csv(
        DATA_DIR / "whole_layout_ripley_profiles.tsv", sep="\t", index=False
    )
    whole_layout_ripley_agreement_table.to_csv(
        DATA_DIR / "whole_layout_ripley_agreement.tsv",
        sep="\t",
        index=False,
    )
    gene_metric_table.to_csv(
        DATA_DIR / "gene_spatial_metrics.tsv", sep="\t", index=False
    )
    gene_agreement_table.to_csv(
        DATA_DIR / "gene_agreement.tsv", sep="\t", index=False
    )
    expression_metric_table.to_csv(
        EXPRESSION_DATA_DIR / "expression_fidelity.tsv",
        sep="\t",
        index=False,
    )
    expression_agreement_table.to_csv(
        EXPRESSION_DATA_DIR / "expression_agreement.tsv",
        sep="\t",
        index=False,
    )
    molecular_replicate_agreement_table.to_csv(
        MOLECULAR_DATA_DIR / "molecular_replicate_agreement.tsv",
        sep="\t",
        index=False,
    )
    colocalization_table.to_csv(
        DATA_DIR / "colocalization_matrices.tsv", sep="\t", index=False
    )
    colocalization_agreement_table.to_csv(
        DATA_DIR / "colocalization_agreement.tsv", sep="\t", index=False
    )
    niche_metric_table.to_csv(
        DATA_DIR / "niche_concordance.tsv", sep="\t", index=False
    )
    niche_composition_table.to_csv(
        DATA_DIR / "niche_composition.tsv", sep="\t", index=False
    )
    summary.to_csv(DATA_DIR / "summary_metrics.tsv", sep="\t", index=False)

    config = {
        "analysis": "R1-2/R1-3 local reference-guided validation",
        "reference_scope": "one 1 mm x 1 mm Xenium breast-tumor tile",
        "fit_reused": True,
        "held_out_data": False,
        "simspace_version": simspace_version,
        "simspace_source_commit": SIMSPACE_SOURCE_COMMIT,
        "simulation_seeds": list(SIMULATION_SEEDS),
        "knn_k": K_NEIGHBORS,
        "ripley_scope": (
            "all cell coordinates in the complete tile, without cell-type "
            "or gene stratification"
        ),
        "ripley_radii_after_axiswise_unit_normalization": (
            WHOLE_LAYOUT_RIPLEY_RADII.tolist()
        ),
        "ripley_edge_correction": "none, matching the Figure 2 implementation",
        "minimum_celltype_count": MIN_CELLTYPE_COUNT,
        "gene_detection_fraction_range": [
            MIN_GENE_DETECTION,
            MAX_GENE_DETECTION,
        ],
        "colocalization_permutations": PERMUTATIONS,
        "gene_diagnostics": (
            "frozen seed-0 molecular realization aligned to the exact "
            "replayed seed-0 coordinates"
        ),
        "expression_fidelity": (
            "Pearson correlation and RMSE between Xenium and the frozen "
            "seed-0 molecular realization for log1p gene-wise raw-count "
            "means and log1p unbiased raw-count variances across all "
            "shared genes"
        ),
        "molecular_replicate_fidelity": (
            "The scDesign3 marginal and copula models were fitted once to "
            "the Xenium tile, then molecular draws for spatial seeds 0--9 "
            "were generated with matching molecular seeds. For each seed, "
            "Pearson correlation and RMSE compare the log1p gene-wise mean "
            "count vector with Xenium across all 220 genes."
        ),
        "molecular_model_fit_seed": 0,
        "molecular_draw_seeds": list(SIMULATION_SEEDS),
        "molecular_model": {
            "marginal_family": "negative binomial",
            "mu_formula": "celltype",
            "sigma_formula": "celltype",
            "copula": "gaussian",
            "correlation_group": "celltype",
            "scDesign3_version": "1.5.0",
        },
        "niche_diagnostics": "frozen Figure 4 BANKSY outputs",
    }
    (DATA_DIR / "analysis_config.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    versions = pd.DataFrame(
        [
            {"software": "python", "version": platform.python_version()},
            {"software": "numpy", "version": np.__version__},
            {"software": "pandas", "version": pd.__version__},
            {"software": "scipy", "version": scipy.__version__},
            {"software": "matplotlib", "version": matplotlib.__version__},
            {"software": "seaborn", "version": sns.__version__},
            {
                "software": "scikit-learn",
                "version": sys.modules["sklearn"].__version__,
            },
            {
                "software": "simspace",
                "version": (
                    f"{simspace_version}; source commit "
                    f"{SIMSPACE_SOURCE_COMMIT}"
                ),
            },
        ]
    )
    versions.to_csv(DATA_DIR / "software_versions.tsv", sep="\t", index=False)

    median_simulated_colocalization = np.median(
        np.stack(simulated_matrices), axis=0
    )
    figure = plot_figure(
        cell_types,
        celltype_agreement_table,
        whole_layout_ripley_table,
        whole_layout_ripley_agreement_table,
        gene_metric_table,
        gene_agreement_table,
        expression_metric_table,
        expression_agreement_table,
        molecular_replicate_agreement_table,
        reference_colocalization,
        median_simulated_colocalization,
        niche_composition_table,
    )
    local_output = SCRIPT_DIR / "SFig9.png"
    example_output = EXAMPLE_DIR / "SFig9.png"
    figure.savefig(local_output, dpi=300, bbox_inches="tight")
    shutil.copyfile(local_output, example_output)
    plt.close(figure)

    print(summary.to_string(index=False))
    print(f"\nWrote {local_output.relative_to(REPO_ROOT)}")
    print(f"Wrote {example_output.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
