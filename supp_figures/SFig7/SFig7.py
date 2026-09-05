"""Supplementary Figure 7: quantitative reference-free 3D simulation.

The analysis uses one fixed SimSpace configuration across ten seeds. It
measures spatial organization directly in three dimensions and evaluates a
small downstream task against exported niche truth. The default entry point
renders from archived results; SFig7_generate_data.py triggers the full
simulation and data-export workflow.
"""

from __future__ import annotations

import json
import os
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/simspace-r1-5-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/simspace-r1-5-cache")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import sklearn
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler

import simspace as ss


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "Panel_A_F_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig7"

SEEDS = tuple(range(10))
SHAPE = (50, 50, 20)
N_NICHES = 3
N_STATES = 6
N_GENES = 36
NICHE_SWEEPS = 4
PHENOTYPE_SWEEPS = 5
PHI = 5.0
DENSITY_RETENTION = 0.85
PERTURBATION_SD = 0.40
MORAN_K = 18
RECOVERY_K = 24
N_INIT = 50
SIX_CONNECTED_NEIGHBORHOOD = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)
SIMSPACE_EXPECTED_VERSION = "0.4.2"
SIMSPACE_SOURCE_COMMIT = "747bb234020f807c8fd9963310cd687dd70f1925"

LEVEL_LABELS = {
    "niche_indicator": "Niche indicators",
    "phenotype_indicator": "Phenotype indicators",
    "marker_gene": "Marker-gene expression",
}
METHOD_LABELS = {
    "expression_only": "Expression only",
    "expression_plus_3d": "Expression + 3D context",
}
PALETTE = ("#2E86AB", "#F18F01", "#6A994E")


def package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "not installed as a distribution"


def fixed_parameter_matrices() -> tuple[np.ndarray, list[np.ndarray]]:
    """Return the pre-specified niche and niche-specific phenotype matrices."""
    niche_theta = np.full((N_NICHES, N_NICHES), 0.05, dtype=float)
    np.fill_diagonal(niche_theta, 1.0)

    theta_list: list[np.ndarray] = []
    for niche in range(N_NICHES):
        theta = np.full((N_STATES, N_STATES), -0.15, dtype=float)
        np.fill_diagonal(theta, 0.25)
        preferred = (2 * niche, 2 * niche + 1)
        theta[np.ix_(preferred, preferred)] = 1.0
        theta_list.append(theta)
    return niche_theta, theta_list


def knn_indices(coordinates: np.ndarray, k: int) -> np.ndarray:
    effective_k = min(k, len(coordinates) - 1)
    neighbors = cKDTree(coordinates).query(
        coordinates, k=effective_k + 1
    )[1]
    if neighbors.ndim == 1:
        neighbors = neighbors[:, None]
    return neighbors[:, 1:]


def moran_i(values: np.ndarray, neighbors: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    denominator = float(centered @ centered)
    if denominator == 0:
        return np.nan
    source = np.repeat(np.arange(len(values)), neighbors.shape[1])
    target = neighbors.reshape(-1)
    return float(
        len(values)
        / len(source)
        * np.sum(centered[source] * centered[target])
        / denominator
    )


def add_exported_niche_truth(simulation: ss.SimSpace) -> None:
    retained = simulation.grid.reshape(-1) >= 0
    niche_truth = simulation.niche.reshape(-1)[retained].astype(int)
    if len(niche_truth) != len(simulation.meta):
        raise RuntimeError("Niche truth is not aligned with exported cells.")
    simulation.meta["niche"] = niche_truth
    simulation.meta["state"] = simulation.meta["state"].astype(int)


def simulate(seed: int) -> ss.SimSpace:
    niche_theta, theta_list = fixed_parameter_matrices()
    simulation = ss.SimSpace(
        shape=SHAPE,
        num_states=N_STATES,
        num_iterations=4,
        theta=theta_list,
        phi=PHI,
        neighborhood=list(SIX_CONNECTED_NEIGHBORHOOD),
        random_seed=seed,
    )
    simulation.initialize3D()
    simulation.create_niche3D(
        num_niches=N_NICHES,
        n_iter=10,
        theta_niche=niche_theta,
        neighborhood=list(SIX_CONNECTED_NEIGHBORHOOD),
    )
    simulation.gibbs_sampler3D()
    simulation.density_sampler(
        np.full(N_NICHES, DENSITY_RETENTION, dtype=float)
    )
    simulation.perturbation3D(step=PERTURBATION_SD)
    add_exported_niche_truth(simulation)
    simulation.create_omics(
        n_genes=N_GENES,
        bg_ratio=1 / 3,
        bg_param=(1.0, 1.0),
        marker_param=(10.0, 0.8),
        lr_ratio=0.0,
        spatial=False,
    )
    return simulation


def spatial_metrics(
    simulation: ss.SimSpace,
    seed: int,
    neighbors: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    niches = simulation.meta["niche"].to_numpy(dtype=int)
    states = simulation.meta["state"].to_numpy(dtype=int)

    for niche in range(N_NICHES):
        rows.append(
            {
                "seed": seed,
                "feature_level": "niche_indicator",
                "component": f"niche_{niche}",
                "moran_i": moran_i((niches == niche).astype(float), neighbors),
            }
        )
    for state in range(N_STATES):
        rows.append(
            {
                "seed": seed,
                "feature_level": "phenotype_indicator",
                "component": f"state_{state}",
                "moran_i": moran_i((states == state).astype(float), neighbors),
            }
        )

    marker_rows = simulation.gene_meta[
        simulation.gene_meta["Marker"].astype(int) >= 0
    ]
    for gene_id in marker_rows["GeneID"].astype(int):
        gene = f"Gene_{gene_id}"
        rows.append(
            {
                "seed": seed,
                "feature_level": "marker_gene",
                "component": gene,
                "moran_i": moran_i(
                    np.log1p(simulation.omics[gene].to_numpy(dtype=float)),
                    neighbors,
                ),
            }
        )
    return pd.DataFrame(rows)


def align_cluster_labels(truth: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    confusion = np.zeros((N_NICHES, N_NICHES), dtype=int)
    for true_label, predicted_label in zip(truth, predicted, strict=True):
        confusion[true_label, predicted_label] += 1
    true_order, predicted_order = linear_sum_assignment(-confusion)
    mapping = dict(zip(predicted_order, true_order, strict=True))
    return np.asarray([mapping[label] for label in predicted], dtype=int)


def recover_niches(
    simulation: ss.SimSpace,
    seed: int,
    neighbors: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    truth = simulation.meta["niche"].to_numpy(dtype=int)
    log_counts = np.log1p(simulation.omics.to_numpy(dtype=float))
    expression = StandardScaler().fit_transform(log_counts)
    local_expression = np.empty_like(log_counts)
    for cell in range(len(log_counts)):
        local_expression[cell] = np.vstack(
            (log_counts[cell], log_counts[neighbors[cell]])
        ).mean(axis=0)
    local_expression = StandardScaler().fit_transform(local_expression)
    expression_plus_3d = np.hstack((expression, local_expression))

    predictions: dict[str, np.ndarray] = {}
    rows: list[dict[str, object]] = []
    for method, features in (
        ("expression_only", expression),
        ("expression_plus_3d", expression_plus_3d),
    ):
        predicted = KMeans(
            n_clusters=N_NICHES,
            n_init=N_INIT,
            random_state=1000 + seed,
        ).fit_predict(features)
        predictions[method] = predicted
        rows.extend(
            [
                {
                    "seed": seed,
                    "method": method,
                    "metric": "ARI",
                    "value": adjusted_rand_score(truth, predicted),
                    "n_cells": len(truth),
                },
                {
                    "seed": seed,
                    "method": method,
                    "metric": "NMI",
                    "value": normalized_mutual_info_score(truth, predicted),
                    "n_cells": len(truth),
                },
            ]
        )
    return (
        pd.DataFrame(rows),
        align_cluster_labels(truth, predictions["expression_only"]),
        align_cluster_labels(truth, predictions["expression_plus_3d"]),
    )


def summarize_spatial_metrics(raw_metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (seed, feature_level), subset in raw_metrics.groupby(
        ["seed", "feature_level"]
    ):
        aggregation = "median" if feature_level == "marker_gene" else "mean"
        value = (
            subset["moran_i"].median()
            if aggregation == "median"
            else subset["moran_i"].mean()
        )
        rows.append(
            {
                "seed": seed,
                "feature_level": feature_level,
                "aggregation": aggregation,
                "n_components": len(subset),
                "moran_i": value,
            }
        )
    return pd.DataFrame(rows)


def centered_offsets(n_values: int, spread: float = 0.075) -> np.ndarray:
    if n_values <= 1:
        return np.zeros(n_values)
    return np.linspace(-spread, spread, n_values)


def plot_recovery_scores(axis: plt.Axes, recovery: pd.DataFrame) -> None:
    method_style = {
        "expression_only": {
            "color": "#8C8C8C",
            "marker": "o",
            "offset": 0.14,
        },
        "expression_plus_3d": {
            "color": "#0072B2",
            "marker": "D",
            "offset": -0.14,
        },
    }
    metric_y = {"ARI": 1.0, "NMI": 0.0}
    for metric, base_y in metric_y.items():
        for method, style in method_style.items():
            values = (
                recovery[
                    (recovery["metric"] == metric)
                    & (recovery["method"] == method)
                ]
                .sort_values("seed")["value"]
                .to_numpy()
            )
            y = base_y + style["offset"]
            axis.hlines(
                y,
                values.min(),
                values.max(),
                color=style["color"],
                linewidth=1.8,
                alpha=0.55,
                zorder=1,
            )
            axis.scatter(
                values,
                y + centered_offsets(len(values)),
                color=style["color"],
                marker=style["marker"],
                s=24,
                alpha=0.78,
                linewidths=0,
                zorder=2,
            )
            axis.scatter(
                np.median(values),
                y,
                color=style["color"],
                marker=style["marker"],
                s=62,
                edgecolor="white",
                linewidth=0.9,
                zorder=3,
            )

    axis.set_yticks([1.0, 0.0], ["ARI", "NMI"])
    axis.set_xlim(0, 1)
    axis.set_ylim(-0.45, 1.45)
    axis.set_xlabel("Recovery score")
    axis.set_ylabel("")
    axis.grid(axis="y", visible=False)
    axis.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker=style["marker"],
                color=style["color"],
                label=METHOD_LABELS[method],
                markerfacecolor=style["color"],
                markersize=5,
                linewidth=1.5,
            )
            for method, style in method_style.items()
        ],
        loc="lower right",
        frameon=False,
        fontsize=7.5,
    )


def plot_recovery_gains(axis: plt.Axes, recovery: pd.DataFrame) -> None:
    metric_style = {
        "ARI": {"color": "#0072B2", "marker": "o", "y": 1.0},
        "NMI": {"color": "#E69F00", "marker": "s", "y": 0.0},
    }
    paired = recovery.pivot(
        index="seed",
        columns=["metric", "method"],
        values="value",
    )
    maximum = 0.0
    for metric, style in metric_style.items():
        gains = (
            paired[(metric, "expression_plus_3d")]
            - paired[(metric, "expression_only")]
        ).to_numpy()
        maximum = max(maximum, float(gains.max()))
        axis.hlines(
            style["y"],
            gains.min(),
            gains.max(),
            color=style["color"],
            linewidth=2.0,
            alpha=0.55,
            zorder=1,
        )
        axis.scatter(
            gains,
            style["y"] + centered_offsets(len(gains), spread=0.11),
            color=style["color"],
            marker=style["marker"],
            s=30,
            alpha=0.82,
            linewidths=0,
            zorder=2,
        )
        median_gain = float(np.median(gains))
        axis.scatter(
            median_gain,
            style["y"],
            color=style["color"],
            marker="D",
            s=68,
            edgecolor="white",
            linewidth=0.9,
            zorder=3,
        )
        axis.text(
            gains.max() + 0.018,
            style["y"],
            rf"median $\Delta$={median_gain:.2f}",
            va="center",
            fontsize=7.5,
            color=style["color"],
        )

    axis.axvline(0, color="#4D4D4D", linewidth=0.9, linestyle="--")
    axis.set_yticks([1.0, 0.0], ["ARI", "NMI"])
    axis.set_xlim(-0.02, max(0.50, maximum + 0.18))
    axis.set_ylim(-0.45, 1.45)
    axis.set_xlabel("Gain from 3D context")
    axis.set_ylabel("")
    axis.grid(axis="y", visible=False)


def plot_figure(
    representative: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    recovery: pd.DataFrame,
) -> plt.Figure:
    np.random.seed(1729)  # Keep seaborn point jitter identical across runs.
    sns.set_theme(style="whitegrid", context="paper")
    figure = plt.figure(figsize=(12.4, 7.3), dpi=300)
    grid = figure.add_gridspec(2, 3, height_ratios=(1.08, 1.0))

    volume_axis = figure.add_subplot(grid[0, 0], projection="3d")
    volume_axis.scatter(
        representative["x"],
        representative["y"],
        representative["z"],
        c=[PALETTE[value] for value in representative["niche"]],
        s=0.75,
        alpha=0.22,
        linewidths=0,
        rasterized=True,
    )
    volume_axis.set_xlabel("x")
    volume_axis.set_ylabel("y")
    volume_axis.set_zlabel("z")
    volume_axis.set_title("A   Exported 3D niche truth", loc="left", fontweight="bold")
    volume_axis.view_init(elev=20, azim=-55)
    volume_axis.set_box_aspect(SHAPE)
    volume_axis.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=PALETTE[niche],
                markeredgecolor="none",
                label=f"Niche {niche + 1}",
                markersize=5,
            )
            for niche in range(N_NICHES)
        ],
        loc="upper left",
        frameon=False,
        fontsize=8,
    )

    central_layer = (SHAPE[2] - 1) / 2
    layer = representative[
        np.abs(representative["z"] - central_layer) < 0.5
    ].copy()
    for column, title, grid_column in (
        ("niche", "B   Truth, central z slab", 1),
        ("expression_plus_3d_cluster", "C   Recovered with 3D context", 2),
    ):
        axis = figure.add_subplot(grid[0, grid_column])
        axis.scatter(
            layer["x"],
            layer["y"],
            c=[PALETTE[value] for value in layer[column].astype(int)],
            s=4.0,
            alpha=0.82,
            linewidths=0,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title, loc="left", fontweight="bold")

    structure_axis = figure.add_subplot(grid[1, 0])
    plot_data = spatial_summary.copy()
    plot_data["feature_label"] = plot_data["feature_level"].map(LEVEL_LABELS)
    order = list(LEVEL_LABELS.values())
    sns.boxplot(
        data=plot_data,
        y="feature_label",
        x="moran_i",
        order=order,
        hue="feature_label",
        palette=dict(zip(order, PALETTE, strict=True)),
        width=0.55,
        fliersize=0,
        legend=False,
        ax=structure_axis,
    )
    sns.stripplot(
        data=plot_data,
        y="feature_label",
        x="moran_i",
        order=order,
        color="#202020",
        size=3.5,
        jitter=0.12,
        ax=structure_axis,
    )
    structure_axis.set_xlabel("3D Moran's I")
    structure_axis.set_ylabel("")
    structure_axis.set_xlim(0, 0.72)
    structure_axis.grid(axis="y", visible=False)
    structure_axis.set_title(
        "D   Volumetric spatial organization",
        loc="left",
        fontweight="bold",
    )

    recovery_axis = figure.add_subplot(grid[1, 1])
    plot_recovery_scores(recovery_axis, recovery)
    recovery_axis.set_title(
        "E   Niche recovery across seeds",
        loc="left",
        fontweight="bold",
    )

    gain_axis = figure.add_subplot(grid[1, 2])
    plot_recovery_gains(gain_axis, recovery)
    gain_axis.set_title(
        "F   Paired improvement",
        loc="left",
        fontweight="bold",
    )

    figure.tight_layout()
    return figure


def write_provenance() -> None:
    simspace_version = package_version("simspace")
    if simspace_version != SIMSPACE_EXPECTED_VERSION:
        raise RuntimeError(
            f"Unexpected SimSpace version: {simspace_version}; "
            f"expected {SIMSPACE_EXPECTED_VERSION}."
        )

    niche_theta, theta_list = fixed_parameter_matrices()
    configuration = {
        "analysis": "R1-5 quantitative reference-free 3D simulation",
        "scope_note": (
            "One fixed reference-free configuration; no "
            "reference-guided 3D fitting or validation."
        ),
        "seeds": list(SEEDS),
        "shape": list(SHAPE),
        "num_niches": N_NICHES,
        "num_states": N_STATES,
        "num_genes": N_GENES,
        "niche_sweeps": NICHE_SWEEPS,
        "phenotype_sweeps": PHENOTYPE_SWEEPS,
        "phi": PHI,
        "neighborhood": [list(offset) for offset in SIX_CONNECTED_NEIGHBORHOOD],
        "density_retention": DENSITY_RETENTION,
        "coordinate_perturbation_sd": PERTURBATION_SD,
        "moran_knn": MORAN_K,
        "recovery_knn": RECOVERY_K,
        "recovery_features": {
            "expression_only": "standardized log1p single-cell expression",
            "expression_plus_3d": (
                "concatenated standardized single-cell expression and "
                "standardized mean log1p expression of each cell plus "
                "its nearest 3D neighbors"
            ),
        },
        "kmeans_clusters": N_NICHES,
        "kmeans_n_init": N_INIT,
        "niche_theta": niche_theta.tolist(),
        "phenotype_theta": [theta.tolist() for theta in theta_list],
        "omics": {
            "background_ratio": 1 / 3,
            "background_gamma_shape_scale": [1.0, 1.0],
            "marker_gamma_shape_scale": [10.0, 0.8],
            "ligand_receptor_ratio": 0.0,
            "direct_spatial_effect": False,
        },
        "simspace_version": simspace_version,
        "simspace_source_commit": SIMSPACE_SOURCE_COMMIT,
    }
    (DATA_DIR / "analysis_config.json").write_text(
        json.dumps(configuration, indent=2) + "\n",
        encoding="utf-8",
    )

    software = pd.DataFrame(
        [
            {"software": "Python", "version": platform.python_version()},
            {"software": "SimSpace", "version": simspace_version},
            {"software": "NumPy", "version": np.__version__},
            {"software": "pandas", "version": pd.__version__},
            {"software": "SciPy", "version": scipy.__version__},
            {"software": "scikit-learn", "version": sklearn.__version__},
            {"software": "matplotlib", "version": matplotlib.__version__},
            {"software": "seaborn", "version": sns.__version__},
            {"software": "platform", "version": platform.platform()},
            {"software": "python_executable", "version": sys.executable},
        ]
    )
    software.to_csv(DATA_DIR / "software_versions.tsv", sep="\t", index=False)


def save_figure(
    representative: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    recovery: pd.DataFrame,
) -> None:
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    figure = plot_figure(representative, spatial_summary, recovery)
    figure.savefig(SCRIPT_DIR / "SFig7.png", dpi=300, bbox_inches="tight")
    figure.savefig(EXAMPLE_DIR / "SFig7.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def load_archived_results() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    paths = {
        "cells": DATA_DIR / "simulated_cells.tsv.gz",
        "spatial": DATA_DIR / "spatial_structure_summary.tsv",
        "recovery": DATA_DIR / "domain_recovery_metrics.tsv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Archived SFig7 results are missing. Run "
            "SFig7_generate_data.py first. Missing: "
            + ", ".join(missing)
        )

    cells = pd.read_csv(paths["cells"], sep="\t")
    representative = cells[cells["seed"] == SEEDS[0]].copy()
    spatial_summary = pd.read_csv(paths["spatial"], sep="\t")
    recovery = pd.read_csv(paths["recovery"], sep="\t")

    required_cell_columns = {
        "seed",
        "x",
        "y",
        "z",
        "niche",
        "expression_plus_3d_cluster",
    }
    missing_columns = required_cell_columns.difference(representative.columns)
    if missing_columns:
        raise ValueError(
            "Archived cell metadata lacks required columns: "
            + ", ".join(sorted(missing_columns))
        )
    expected_seeds = set(SEEDS)
    if set(spatial_summary["seed"]) != expected_seeds:
        raise ValueError("Archived spatial summaries do not contain seeds 0--9.")
    if set(recovery["seed"]) != expected_seeds:
        raise ValueError("Archived recovery metrics do not contain seeds 0--9.")
    return representative, spatial_summary, recovery


def generate_data() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    write_provenance()

    spatial_frames: list[pd.DataFrame] = []
    recovery_frames: list[pd.DataFrame] = []
    cell_frames: list[pd.DataFrame] = []
    count_frames: list[pd.DataFrame] = []

    for seed in SEEDS:
        print(f"seed={seed}", flush=True)
        simulation = simulate(seed)
        coordinates = simulation.meta[["x", "y", "z"]].to_numpy(dtype=float)
        moran_neighbors = knn_indices(coordinates, MORAN_K)
        recovery_neighbors = knn_indices(coordinates, RECOVERY_K)
        spatial_frames.append(
            spatial_metrics(simulation, seed, moran_neighbors)
        )
        recovery, expression_cluster, expression_plus_3d_cluster = recover_niches(
            simulation, seed, recovery_neighbors
        )
        recovery_frames.append(recovery)

        cells = simulation.meta[["x", "y", "z", "state", "niche"]].copy()
        cells.insert(0, "cell_id", np.arange(len(cells)))
        cells.insert(0, "seed", seed)
        cells["expression_cluster"] = expression_cluster
        cells["expression_plus_3d_cluster"] = expression_plus_3d_cluster
        cell_frames.append(cells)

        counts = simulation.omics.copy()
        counts.insert(0, "cell_id", np.arange(len(counts)))
        counts.insert(0, "seed", seed)
        count_frames.append(counts)

    raw_spatial = pd.concat(spatial_frames, ignore_index=True)
    spatial_summary = summarize_spatial_metrics(raw_spatial)
    recovery = pd.concat(recovery_frames, ignore_index=True)
    cells = pd.concat(cell_frames, ignore_index=True)
    counts = pd.concat(count_frames, ignore_index=True)

    raw_spatial.to_csv(
        DATA_DIR / "spatial_structure_metrics.tsv", sep="\t", index=False
    )
    spatial_summary.to_csv(
        DATA_DIR / "spatial_structure_summary.tsv", sep="\t", index=False
    )
    recovery.to_csv(
        DATA_DIR / "domain_recovery_metrics.tsv", sep="\t", index=False
    )
    summary_metrics = pd.concat(
        [
            spatial_summary.rename(
                columns={
                    "feature_level": "group",
                    "moran_i": "value",
                }
            ).assign(
                analysis="spatial_structure",
                metric="Moran_I",
            )[
                [
                    "seed",
                    "analysis",
                    "metric",
                    "group",
                    "value",
                    "n_components",
                    "aggregation",
                ]
            ],
            recovery.rename(
                columns={
                    "method": "group",
                    "n_cells": "n_components",
                }
            ).assign(
                analysis="domain_recovery",
                aggregation="none",
            )[
                [
                    "seed",
                    "analysis",
                    "metric",
                    "group",
                    "value",
                    "n_components",
                    "aggregation",
                ]
            ],
        ],
        ignore_index=True,
    )
    summary_metrics.to_csv(
        DATA_DIR / "summary_metrics.tsv", sep="\t", index=False
    )
    cells.to_csv(
        DATA_DIR / "simulated_cells.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    counts.to_csv(
        DATA_DIR / "simulated_counts.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    representative, archived_spatial, archived_recovery = load_archived_results()
    save_figure(representative, archived_spatial, archived_recovery)

    print("\nSpatial structure summary:")
    print(
        spatial_summary.groupby("feature_level")["moran_i"]
        .agg(["median", "min", "max"])
        .to_string()
    )
    print("\nDomain recovery summary:")
    print(
        recovery.groupby(["metric", "method"])["value"]
        .agg(["median", "min", "max"])
        .to_string()
    )


def main() -> None:
    representative, spatial_summary, recovery = load_archived_results()
    save_figure(representative, spatial_summary, recovery)
    print(
        "Rendered SFig7 from archived ten-seed results. "
        "Run SFig7_generate_data.py to regenerate the simulations."
    )


if __name__ == "__main__":
    main()
