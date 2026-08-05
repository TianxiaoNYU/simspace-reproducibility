#!/usr/bin/env python3
"""Render the provisional new Figure 5 spatial-domain benchmark panels."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "simspace-domain-figure-mpl")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "simspace-domain-figure-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
import seaborn as sns


SCRIPT_DIR = Path(__file__).resolve().parent
SUMMARY_PATH = SCRIPT_DIR / "Panel_C_data" / "experiment_summary.tsv"
MAP_CACHE_PATH = SCRIPT_DIR / "Panel_A_B_data" / "representative_maps.tsv.gz"
RAW_DATA_ROOT = SCRIPT_DIR / "Panel_A_B_C_src" / "experiment_data"
RAW_RESULTS_ROOT = SCRIPT_DIR / "Panel_A_B_C_src" / "experiment_results"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR

REPRESENTATIVE_SEED = 0
REPRESENTATIVE_DIFFICULTY = "moderate"
TOPOLOGY_ORDER = ("curved_layers", "irregular_mrf")
TOPOLOGY_LABELS = {
    "curved_layers": "Curved layers",
    "irregular_mrf": "Irregular MRF",
}
DIFFICULTY_ORDER = ("hard", "moderate", "easy")
FOLD_CHANGES = {"hard": 2.0, "moderate": 3.0, "easy": 4.0}

METHOD_ORDER = (
    "graphst",
    "stagate",
    "banksy",
    "spagcn",
    "spclue",
    "stlearn",
    "expression_pca_kmeans",
)
METHOD_LABELS = {
    "graphst": "GraphST",
    "stagate": "STAGATE",
    "banksy": "BANKSY",
    "spagcn": "SpaGCN",
    "spclue": "spCLUE",
    "stlearn": "stLearn",
    "expression_pca_kmeans": "Expression-only",
}

TAB10 = sns.color_palette("tab10", 10)
METHOD_COLORS = {
    "graphst": TAB10[0],
    "stagate": TAB10[1],
    "banksy": TAB10[2],
    "spagcn": TAB10[3],
    "spclue": TAB10[4],
    "stlearn": TAB10[5],
    "expression_pca_kmeans": (0.45, 0.45, 0.45),
}
DOMAIN_COLORS = tuple(TAB10[index] for index in range(4))
METHOD_MARKERS = {
    "graphst": "o",
    "stagate": "X",
    "banksy": "s",
    "spagcn": "D",
    "spclue": "v",
    "stlearn": "P",
    "expression_pca_kmeans": "^",
}
METHOD_LINESTYLES = {
    "graphst": "-",
    "stagate": "--",
    "banksy": ":",
    "spagcn": "-.",
    "spclue": "--",
    "stlearn": ":",
    "expression_pca_kmeans": "-.",
}

PANEL_B_COLUMNS = (
    "truth",
    "expression_pca_kmeans",
    "graphst",
    "stagate",
    "banksy",
    "spagcn",
    "spclue",
    "stlearn",
)


def configure_style() -> None:
    sns.set_style("white")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.5,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.linewidth": 0.9,
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "text.color": "black",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_summary() -> pd.DataFrame:
    summary = pd.read_csv(SUMMARY_PATH, sep="\t")
    required = {
        "dataset",
        "topology",
        "difficulty",
        "domain_fold_change",
        "layout_seed",
        "method",
        "status",
        "adjusted_rand_index",
        "boundary_f1",
    }
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"Summary is missing columns: {sorted(missing)}")
    if len(summary) != 126:
        raise ValueError(f"Expected 126 benchmark rows, found {len(summary)}")
    return summary


def _read_representative_dataset(topology: str) -> pd.DataFrame:
    dataset = f"{topology}_seed{REPRESENTATIVE_SEED}_{REPRESENTATIVE_DIFFICULTY}"
    dataset_dir = RAW_DATA_ROOT / dataset
    result_dir = RAW_RESULTS_ROOT / dataset

    coordinates = pd.read_csv(dataset_dir / "coordinates.tsv", sep="\t")
    truth = pd.read_csv(dataset_dir / "truth.tsv", sep="\t")[
        ["location_id", "domain_truth", "cell_type_truth"]
    ]
    merged = coordinates.merge(truth, on="location_id", validate="one_to_one")
    if len(merged) != 3000:
        raise ValueError(f"Expected 3,000 locations in {dataset}, found {len(merged)}")

    for method in METHOD_ORDER:
        assignments = pd.read_csv(
            result_dir / method / "assignments.tsv", sep="\t"
        ).rename(columns={"predicted_domain": method})
        merged = merged.merge(
            assignments[["location_id", method]],
            on="location_id",
            validate="one_to_one",
        )
    merged.insert(0, "dataset", dataset)
    merged.insert(1, "topology", topology)
    return merged


def build_representative_cache(force: bool = False) -> pd.DataFrame:
    if MAP_CACHE_PATH.exists() and not force:
        return pd.read_csv(MAP_CACHE_PATH, sep="\t")
    frames = [_read_representative_dataset(topology) for topology in TOPOLOGY_ORDER]
    cache = pd.concat(frames, ignore_index=True)
    if len(cache) != 6000:
        raise ValueError(f"Expected 6,000 cached locations, found {len(cache)}")
    MAP_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    cache.to_csv(
        MAP_CACHE_PATH,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.8g",
    )
    return cache


def load_representative_cache(force_rebuild: bool = False) -> pd.DataFrame:
    cache = build_representative_cache(force=force_rebuild)
    required = {
        "dataset",
        "topology",
        "location_id",
        "x",
        "y",
        "domain_truth",
        "cell_type_truth",
        *METHOD_ORDER,
    }
    missing = required.difference(cache.columns)
    if missing:
        raise ValueError(f"Map cache is missing columns: {sorted(missing)}")
    if len(cache) != 6000 or cache.isna().any().any():
        raise ValueError("Representative-map cache is incomplete")
    return cache


def hungarian_align(predicted: pd.Series, truth: pd.Series) -> np.ndarray:
    pred_values = predicted.astype(str).to_numpy()
    truth_values = truth.astype(int).to_numpy()
    pred_labels = sorted(pd.unique(pred_values).tolist())
    truth_labels = sorted(pd.unique(truth_values).tolist())
    contingency = np.zeros((len(pred_labels), len(truth_labels)), dtype=int)
    for row, pred_label in enumerate(pred_labels):
        selected = pred_values == pred_label
        for col, truth_label in enumerate(truth_labels):
            contingency[row, col] = int(np.sum(selected & (truth_values == truth_label)))
    rows, cols = linear_sum_assignment(-contingency)
    mapping = {pred_labels[row]: truth_labels[col] for row, col in zip(rows, cols)}
    unused = [label for label in truth_labels if label not in mapping.values()]
    for pred_label in pred_labels:
        if pred_label not in mapping:
            mapping[pred_label] = unused.pop(0) if unused else truth_labels[0]
    return np.asarray([mapping[value] for value in pred_values], dtype=int)


def scatter_labels(
    ax: plt.Axes,
    frame: pd.DataFrame,
    labels: np.ndarray,
    colors: tuple | list,
    point_size: float,
) -> None:
    for label_index, color in enumerate(colors):
        selected = labels == label_index
        ax.scatter(
            frame.loc[selected, "x"],
            frame.loc[selected, "y"],
            s=point_size,
            color=color,
            edgecolors="none",
            rasterized=True,
        )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def save_panel(fig: plt.Figure, output_dir: Path, stem: str, dpi: int) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{stem}.png"
    svg_path = output_dir / f"{stem}.svg"
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(
        svg_path,
        bbox_inches="tight",
        pad_inches=0.03,
        metadata={"Date": None, "Creator": "Fig5_new.py"},
    )
    plt.close(fig)
    return [png_path, svg_path]


def build_composition_matrix(cache: pd.DataFrame) -> pd.DataFrame:
    composition = (
        cache.groupby(["topology", "domain_truth", "cell_type_truth"])
        .size()
        .rename("count")
        .reset_index()
    )
    composition["proportion"] = composition["count"] / composition.groupby(
        ["topology", "domain_truth"]
    )["count"].transform("sum")
    return (
        composition.groupby(["domain_truth", "cell_type_truth"])["proportion"]
        .mean()
        .unstack(fill_value=0.0)
        .reindex(index=range(4), columns=range(10), fill_value=0.0)
    )


def draw_panel_a(
    container: plt.Figure,
    cache: pd.DataFrame,
    panel_label: str | None = None,
) -> None:
    """Draw two reference-free geometries above their composition heatmap."""

    grid = container.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.92,
        top=0.88,
        bottom=0.10,
        height_ratios=(1.0, 0.70),
        wspace=0.12,
        hspace=0.72,
    )
    map_axes = [container.add_subplot(grid[0, index]) for index in range(2)]
    heat_ax = container.add_subplot(grid[1, :])

    for ax, topology in zip(map_axes, TOPOLOGY_ORDER):
        frame = cache.loc[cache["topology"] == topology]
        scatter_labels(
            ax,
            frame,
            frame["domain_truth"].astype(int).to_numpy(),
            DOMAIN_COLORS,
            point_size=2.5,
        )
        ax.set_anchor("N")

    # Fixed container coordinates keep the titles aligned even though the two
    # equal-aspect maps have different spatial extents.
    container.text(0.275, 0.925, TOPOLOGY_LABELS[TOPOLOGY_ORDER[0]], ha="center")
    container.text(0.725, 0.925, TOPOLOGY_LABELS[TOPOLOGY_ORDER[1]], ha="center")

    matrix = build_composition_matrix(cache)
    image = heat_ax.imshow(
        matrix.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=0.15,
    )
    heat_ax.set_title("Cell-type composition")
    heat_ax.set_xlabel("Cell type")
    heat_ax.set_ylabel("Spatial domain")
    heat_ax.set_xticks(np.arange(10), np.arange(1, 11))
    heat_ax.set_yticks(np.arange(4), np.arange(1, 5))
    colorbar = container.colorbar(image, ax=heat_ax, fraction=0.030, pad=0.025)
    colorbar.set_label("Proportion")

    domain_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=5.5,
            markerfacecolor=color,
            markeredgecolor="none",
            label=f"Domain {index + 1}",
        )
        for index, color in enumerate(DOMAIN_COLORS)
    ]
    container.legend(
        handles=domain_handles,
        loc="center",
        bbox_to_anchor=(0.50, 0.505),
        ncol=4,
        frameon=False,
        handletextpad=0.3,
        columnspacing=1.0,
    )
    if panel_label is not None:
        container.text(
            0.005,
            0.995,
            panel_label,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
        )


def make_panel_a(cache: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(6.2, 5.65))
    draw_panel_a(fig, cache)
    return fig


def draw_panel_b(
    container: plt.Figure,
    cache: pd.DataFrame,
    summary: pd.DataFrame,
    panel_label: str | None = None,
) -> None:
    """Representative truth and prediction maps with per-map ARI labels."""

    grid = container.add_gridspec(
        2,
        len(PANEL_B_COLUMNS),
        left=0.075,
        right=0.99,
        top=0.84,
        bottom=0.22,
        wspace=0.06,
        hspace=0.34,
    )
    axes = np.asarray(
        [
            [
                container.add_subplot(grid[row, col])
                for col in range(len(PANEL_B_COLUMNS))
            ]
            for row in range(2)
        ]
    )

    for row, topology in enumerate(TOPOLOGY_ORDER):
        frame = cache.loc[cache["topology"] == topology]
        truth = frame["domain_truth"].astype(int)
        for col, item in enumerate(PANEL_B_COLUMNS):
            ax = axes[row, col]
            if item == "truth":
                labels = truth.to_numpy()
                title = "Truth"
            else:
                labels = hungarian_align(frame[item], truth)
                title = METHOD_LABELS[item]
            scatter_labels(ax, frame, labels, DOMAIN_COLORS, point_size=2.5)
            if item == "truth":
                ari = 1.0
            else:
                dataset = (
                    f"{topology}_seed{REPRESENTATIVE_SEED}_"
                    f"{REPRESENTATIVE_DIFFICULTY}"
                )
                result = summary.loc[
                    (summary["dataset"] == dataset) & (summary["method"] == item),
                    ["status", "adjusted_rand_index"],
                ]
                if len(result) != 1 or result.iloc[0]["status"] != "success":
                    raise ValueError(
                        f"Missing valid representative ARI for {dataset}/{item}"
                    )
                ari = float(result.iloc[0]["adjusted_rand_index"])
            ax.text(
                0.5,
                -0.10,
                f"ARI: {ari:.2f}",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=8.0,
            )
            if row == 0:
                ax.set_title(title, pad=4.0)
            if col == 0:
                ax.set_ylabel(TOPOLOGY_LABELS[topology], fontweight="normal")

    domain_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=5.5,
            markerfacecolor=color,
            markeredgecolor="none",
            label=f"Domain {index + 1}",
        )
        for index, color in enumerate(DOMAIN_COLORS)
    ]
    container.legend(
        handles=domain_handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.005),
        ncol=4,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.5,
    )
    if panel_label is not None:
        container.text(
            0.005,
            0.995,
            panel_label,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
        )


def make_panel_b(cache: pd.DataFrame, summary: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(12.5, 4.0))
    draw_panel_b(fig, cache, summary)
    return fig


def _method_curve(
    summary: pd.DataFrame,
    topology: str,
    method: str,
    metric: str,
) -> pd.DataFrame:
    records: list[dict[str, float]] = []
    for difficulty in DIFFICULTY_ORDER:
        values = summary.loc[
            (summary["topology"] == topology)
            & (summary["method"] == method)
            & (summary["difficulty"] == difficulty)
            & (summary["status"] == "success"),
            metric,
        ].dropna()
        if values.empty:
            continue
        records.append(
            {
                "fold_change": FOLD_CHANGES[difficulty],
                "median": float(values.median()),
                "minimum": float(values.min()),
                "maximum": float(values.max()),
            }
        )
    return pd.DataFrame.from_records(records)


def _draw_metric_curves(
    ax: plt.Axes,
    summary: pd.DataFrame,
    topology: str,
    metric: str,
    y_limits: tuple[float, float],
) -> None:
    for method in METHOD_ORDER:
        curve = _method_curve(summary, topology, method, metric)
        x = curve["fold_change"].to_numpy()
        median = curve["median"].to_numpy()
        ax.fill_between(
            x,
            curve["minimum"].to_numpy(),
            curve["maximum"].to_numpy(),
            color=METHOD_COLORS[method],
            alpha=0.10,
            linewidth=0,
        )
        ax.plot(
            x,
            median,
            color=METHOD_COLORS[method],
            marker=METHOD_MARKERS[method],
            linestyle=METHOD_LINESTYLES[method],
            linewidth=1.8,
            markersize=5.2,
            markeredgecolor="white",
            markeredgewidth=0.55,
            label=METHOD_LABELS[method],
        )
    ax.set_xlim(1.85, 4.15)
    ax.set_ylim(*y_limits)
    ax.set_xticks([2, 3, 4], ["2×\nHard", "3×\nModerate", "4×\nEasy"])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])


def draw_panel_c(
    container: plt.Figure,
    summary: pd.DataFrame,
    panel_label: str | None = None,
) -> None:
    """Combined ARI and boundary-F1 curves with one shared legend."""

    grid = container.add_gridspec(
        2,
        2,
        left=0.08,
        right=0.78,
        top=0.92,
        bottom=0.14,
        wspace=0.18,
        hspace=0.16,
    )
    axes = np.empty((2, 2), dtype=object)
    for row in range(2):
        for col in range(2):
            axes[row, col] = container.add_subplot(
                grid[row, col],
                sharex=axes[0, col] if row == 1 else None,
                sharey=axes[row, 0] if col == 1 else None,
            )

    metrics = (
        ("adjusted_rand_index", "Adjusted Rand index", (-0.055, 1.025)),
        ("boundary_f1", "Boundary F1", (-0.025, 1.025)),
    )
    for row, (metric, y_label, y_limits) in enumerate(metrics):
        for col, topology in enumerate(TOPOLOGY_ORDER):
            ax = axes[row, col]
            _draw_metric_curves(ax, summary, topology, metric, y_limits)
            if row == 0:
                ax.set_title(TOPOLOGY_LABELS[topology])
                ax.tick_params(axis="x", labelbottom=False)
            if col == 1:
                ax.tick_params(axis="y", left=False, labelleft=False)
        axes[row, 0].set_ylabel(y_label)
    container.supxlabel("Domain-specific expression strength", x=0.43, y=0.025)

    handles, labels = axes[0, 1].get_legend_handles_labels()
    container.legend(
        handles,
        labels,
        title="Method",
        loc="center left",
        bbox_to_anchor=(0.80, 0.52),
        frameon=True,
    )
    if panel_label is not None:
        container.text(
            0.005,
            0.995,
            panel_label,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
        )


def make_panel_c(summary: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=(8.8, 5.65))
    draw_panel_c(fig, summary)
    return fig


def make_complete_figure(cache: pd.DataFrame, summary: pd.DataFrame) -> plt.Figure:
    """Assemble A and C above the full-width representative-map panel B."""

    fig = plt.figure(figsize=(12.5, 9.45), facecolor="white")
    outer = fig.add_gridspec(
        2,
        1,
        left=0.045,
        right=0.99,
        top=0.965,
        bottom=0.045,
        height_ratios=(1.25, 1.0),
        hspace=0.20,
    )
    top = outer[0].subgridspec(
        1,
        2,
        width_ratios=(0.43, 0.57),
        wspace=0.24,
    )

    # Panel A: fixed title row, spatial maps, legend, then the full-width heatmap.
    panel_a_grid = top[0, 0].subgridspec(
        4,
        3,
        height_ratios=(0.12, 1.0, 0.16, 0.62),
        width_ratios=(1.0, 1.0, 0.065),
        wspace=0.12,
        hspace=0.14,
    )
    for col, topology in enumerate(TOPOLOGY_ORDER):
        title_ax = fig.add_subplot(panel_a_grid[0, col])
        title_ax.axis("off")
        title_ax.text(
            0.5,
            0.35,
            TOPOLOGY_LABELS[topology],
            ha="center",
            va="center",
            fontsize=10.5,
        )
        map_ax = fig.add_subplot(panel_a_grid[1, col])
        frame = cache.loc[cache["topology"] == topology]
        scatter_labels(
            map_ax,
            frame,
            frame["domain_truth"].astype(int).to_numpy(),
            DOMAIN_COLORS,
            point_size=2.5,
        )
        map_ax.set_anchor("N")

    domain_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=5.5,
            markerfacecolor=color,
            markeredgecolor="none",
            label=f"Domain {index + 1}",
        )
        for index, color in enumerate(DOMAIN_COLORS)
    ]
    panel_a_legend = fig.add_subplot(panel_a_grid[2, :2])
    panel_a_legend.axis("off")
    panel_a_legend.legend(
        handles=domain_handles,
        loc="center",
        ncol=4,
        frameon=False,
        handletextpad=0.3,
        columnspacing=1.0,
    )

    heat_ax = fig.add_subplot(panel_a_grid[3, :2])
    heat_image = heat_ax.imshow(
        build_composition_matrix(cache).to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=0.0,
        vmax=0.20,
    )
    heat_ax.set_title("Cell-type composition")
    heat_ax.set_xlabel("Cell type")
    heat_ax.set_ylabel("Spatial domain")
    heat_ax.set_xticks(np.arange(10), np.arange(1, 11))
    heat_ax.set_yticks(np.arange(4), np.arange(1, 5))
    heat_colorbar_ax = fig.add_subplot(panel_a_grid[3, 2])
    heat_colorbar = fig.colorbar(heat_image, cax=heat_colorbar_ax)
    heat_colorbar.ax.set_title("Proportion", fontsize=8.5, pad=4.0)

    # Panel C: the 2 x 2 performance grid, a dedicated legend column, and one
    # shared signal-strength label row.
    panel_c_grid = top[0, 1].subgridspec(
        3,
        3,
        height_ratios=(1.0, 1.0, 0.12),
        width_ratios=(1.0, 1.0, 0.54),
        wspace=0.22,
        hspace=0.18,
    )
    metric_axes = np.empty((2, 2), dtype=object)
    metrics = (
        ("adjusted_rand_index", "Adjusted Rand index", (-0.055, 1.025)),
        ("boundary_f1", "Boundary F1", (-0.025, 1.025)),
    )
    for row, (metric, y_label, y_limits) in enumerate(metrics):
        for col, topology in enumerate(TOPOLOGY_ORDER):
            metric_axes[row, col] = fig.add_subplot(
                panel_c_grid[row, col],
                sharex=metric_axes[0, col] if row == 1 else None,
                sharey=metric_axes[row, 0] if col == 1 else None,
            )
            ax = metric_axes[row, col]
            _draw_metric_curves(ax, summary, topology, metric, y_limits)
            if row == 0:
                ax.set_title(TOPOLOGY_LABELS[topology])
                ax.tick_params(axis="x", labelbottom=False)
            if col == 1:
                ax.tick_params(axis="y", left=False, labelleft=False)
        metric_axes[row, 0].set_ylabel(y_label)

    method_handles, method_labels = metric_axes[0, 1].get_legend_handles_labels()
    method_legend_ax = fig.add_subplot(panel_c_grid[:2, 2])
    method_legend_ax.axis("off")
    method_legend_ax.legend(
        method_handles,
        method_labels,
        title="Method",
        loc="center",
        frameon=True,
    )
    signal_label_ax = fig.add_subplot(panel_c_grid[2, :2])
    signal_label_ax.axis("off")
    signal_label_ax.text(
        0.5,
        0.08,
        "Domain-specific expression strength",
        ha="center",
        va="center",
        fontsize=10.0,
    )

    # Panel B: the wide representative-map grid spans the complete second row.
    panel_b_grid = outer[1, 0].subgridspec(
        3,
        len(PANEL_B_COLUMNS),
        height_ratios=(1.0, 1.0, 0.17),
        wspace=0.06,
        hspace=0.34,
    )
    panel_b_axes = np.empty((2, len(PANEL_B_COLUMNS)), dtype=object)
    for row, topology in enumerate(TOPOLOGY_ORDER):
        frame = cache.loc[cache["topology"] == topology]
        truth = frame["domain_truth"].astype(int)
        for col, item in enumerate(PANEL_B_COLUMNS):
            ax = fig.add_subplot(panel_b_grid[row, col])
            panel_b_axes[row, col] = ax
            if item == "truth":
                labels = truth.to_numpy()
                title = "Truth"
                ari = 1.0
            else:
                labels = hungarian_align(frame[item], truth)
                title = METHOD_LABELS[item]
                dataset = (
                    f"{topology}_seed{REPRESENTATIVE_SEED}_"
                    f"{REPRESENTATIVE_DIFFICULTY}"
                )
                result = summary.loc[
                    (summary["dataset"] == dataset) & (summary["method"] == item),
                    ["status", "adjusted_rand_index"],
                ]
                if len(result) != 1 or result.iloc[0]["status"] != "success":
                    raise ValueError(
                        f"Missing valid representative ARI for {dataset}/{item}"
                    )
                ari = float(result.iloc[0]["adjusted_rand_index"])
            scatter_labels(ax, frame, labels, DOMAIN_COLORS, point_size=2.5)
            ax.text(
                0.5,
                -0.10,
                f"ARI: {ari:.2f}",
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=8.0,
            )
            if row == 0:
                ax.set_title(title, pad=4.0)
            if col == 0:
                ax.set_ylabel(TOPOLOGY_LABELS[topology], fontweight="normal")

    panel_b_legend = fig.add_subplot(panel_b_grid[2, :])
    panel_b_legend.axis("off")
    panel_b_legend.legend(
        handles=domain_handles,
        loc="center",
        ncol=4,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.5,
    )

    panel_a_box = top[0, 0].get_position(fig)
    panel_c_box = top[0, 1].get_position(fig)
    panel_b_box = outer[1, 0].get_position(fig)
    label_offsets = {"A": 0.035, "C": 0.055, "B": 0.035}
    for label, box in (("A", panel_a_box), ("C", panel_c_box), ("B", panel_b_box)):
        fig.text(
            box.x0 - label_offsets[label],
            box.y1 + 0.012,
            label,
            ha="left",
            va="top",
            fontsize=14,
            fontweight="bold",
        )
    return fig


def write_manifest(output_dir: Path, outputs: dict[str, list[Path]]) -> Path:
    manifest = {
        "schema_version": "4.0",
        "representative_condition": {
            "difficulty": REPRESENTATIVE_DIFFICULTY,
            "layout_seed": REPRESENTATIVE_SEED,
            "visual_cluster_alignment": "maximum-overlap Hungarian assignment; display only",
        },
        "curve_summary": "median with observed seed range",
        "boundary_summary": "median with observed seed range",
        "combined_performance_panel": {
            "columns": list(TOPOLOGY_ORDER),
            "rows": ["adjusted_rand_index", "boundary_f1"],
            "shared_legend": True,
        },
        "complete_figure_layout": {
            "top_row": ["A", "C"],
            "bottom_row": ["B"],
            "panel_b_spans_full_width": True,
        },
        "signal_axis": "domain-specific expression strength (2x, 3x, 4x)",
        "composition_color_scale": [0.0, 0.20],
        "map_annotation": "ARI for the displayed seed and difficulty",
        "validity_reporting": "supporting documentation only; not drawn in panels",
        "independent_unit": "layout seed",
        "source_files": {
            "experiment_summary.tsv": sha256(SUMMARY_PATH),
            "representative_maps.tsv.gz": sha256(MAP_CACHE_PATH),
        },
        "outputs": {
            panel: {path.name: sha256(path) for path in paths}
            for panel, paths in sorted(outputs.items())
        },
    }
    manifest_path = output_dir / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panels",
        nargs="+",
        choices=("A", "B", "C"),
        default=("A", "B", "C"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--rebuild-map-cache", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()
    requested = tuple(dict.fromkeys(args.panels))
    summary = load_summary()
    cache = load_representative_cache(args.rebuild_map_cache)
    outputs: dict[str, list[Path]] = {}

    if "A" in requested:
        outputs["A"] = save_panel(
            make_panel_a(cache), args.output_dir, "Panel_A", args.dpi
        )
    if "B" in requested:
        outputs["B"] = save_panel(
            make_panel_b(cache, summary),
            args.output_dir,
            "Panel_B",
            args.dpi,
        )
    if "C" in requested:
        outputs["C"] = save_panel(
            make_panel_c(summary), args.output_dir, "Panel_C", args.dpi
        )
    if set(requested) == {"A", "B", "C"}:
        outputs["complete"] = save_panel(
            make_complete_figure(cache, summary),
            args.output_dir,
            "Fig5_new",
            args.dpi,
        )
    manifest_path = write_manifest(args.output_dir, outputs)
    rendered = [*requested]
    if "complete" in outputs:
        rendered.append("complete figure")
    print(f"Rendered {', '.join(rendered)} to {args.output_dir}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
