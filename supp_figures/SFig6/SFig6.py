#!/usr/bin/env python3
"""Supplementary Figure 6: reference-free mouse-cortex comparison.

The analysis compares the fixed reference-free Figure 3 cortex simulation
with a compact mouse visual-cortex STARmap field. STARmap is used only after
simulation and is not used to fit or tune the simulated anchor.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/simspace-r2-4-starmap-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/simspace-r2-4-starmap-cache")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import pearsonr, spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "Panel_A_E_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig6"
RAW_ARCHIVE = DATA_DIR / "STARmap_mouse_visual_cortex.zip"
FIG3_DIR = REPO_ROOT / "main_figures" / "Fig3"
SIMULATED_META = FIG3_DIR / "Panel_A_C_data" / "anchor_cells.tsv.gz"
SIMULATED_EXPRESSION = FIG3_DIR / "Panel_C_data" / "simulated_expression.tsv.gz"

SOURCE_RECORD = "https://zenodo.org/records/10698912"
SOURCE_URL = (
    "https://zenodo.org/records/10698912/files/"
    "STARmap_mouse_visual_cortex.zip?download=1"
)
SOURCE_LICENSE = "CC BY 4.0"
SOURCE_DOI = "10.5281/zenodo.10698912"
ORIGINAL_STARMAP_DOI = "10.1126/science.aat5691"
EXPECTED_RAW_SHA256 = "c1ffc459cf149b25773411fc6c5ef32845104dd8974802f627171fbf394e083c"
EXPECTED_RAW_MD5 = "8f8464ab54504d8854a9fb39cdbcb478"
EXPECTED_RAW_BYTES = 471_002
H5AD_MEMBER = "STARmap_mouse_visual_cortex/STARmap_20180505_BY3_1k.h5ad"
ANNOTATION_MEMBER = (
    "STARmap_mouse_visual_cortex/gt/Annotation_STARmap_20180505_BY3_1k.txt"
)

CLASS_ORDER = ("superficial", "middle", "deep", "inner")
CLASS_LABELS = {
    "superficial": "Superficial",
    "middle": "Middle",
    "deep": "Deep",
    "inner": "Inner",
}
CLASS_COLORS = {
    "superficial": "#56B4E9",
    "middle": "#0072B2",
    "deep": "#D55E00",
    "inner": "#8C6BB1",
}
MARKERS = ("Cux2", "Rorb", "Fezf2", "Deptor", "Tshz2", "Foxp2", "Sulf2", "Ctgf")
REPRESENTATIVE_MARKERS = ("Cux2", "Rorb")
MARKER_EXPECTED_CLASS = {
    "Cux2": "superficial",
    "Rorb": "middle",
    "Fezf2": "deep",
    "Deptor": "deep",
    "Tshz2": "inner",
    "Foxp2": "inner",
    "Sulf2": "inner",
    "Ctgf": "inner",
}

OBSERVED_CELL_CLASS = {
    "eL2/3": "superficial",
    "eL4": "middle",
    "eL5": "deep",
    "eL6-1": "inner",
    "eL6-2": "inner",
}
OBSERVED_DOMAIN = {
    "L2/3": "superficial",
    "L4": "middle",
    "L5": "deep",
    "L6": "inner",
}
SIMULATED_CELL_CLASS = {
    "L2/3 IT": "superficial",
    "L4/5 IT": "middle",
    "L5 ET": "deep",
    "L5 IT": "deep",
    "L5/6 NP": "inner",
    "L6 CT": "inner",
    "L6 IT": "inner",
    "L6 IT Car3": "inner",
    "L6b": "inner",
}
SIMULATED_DOMAIN = {
    "superficial (L2/3)": "superficial",
    "middle (L4/5)": "middle",
    "deep (L5)": "deep",
    "inner (L5/6-L6b)": "inner",
}

EXPRESSION_CMAP = LinearSegmentedColormap.from_list(
    "starmap_expression",
    ["#ECEFF1", "#B8D8E8", "#56B1C9", "#147D92", "#073B4C"],
)


def checksum(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values]
    )


def verify_inputs() -> None:
    missing = [
        str(path)
        for path in (RAW_ARCHIVE, SIMULATED_META, SIMULATED_EXPRESSION)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Missing inputs: "
            + ", ".join(missing)
            + ". Run download_starmap.py and Fig3.py as needed."
        )
    observed = {
        "bytes": RAW_ARCHIVE.stat().st_size,
        "sha256": checksum(RAW_ARCHIVE),
        "md5": checksum(RAW_ARCHIVE, "md5"),
    }
    expected = {
        "bytes": EXPECTED_RAW_BYTES,
        "sha256": EXPECTED_RAW_SHA256,
        "md5": EXPECTED_RAW_MD5,
    }
    if observed != expected:
        raise RuntimeError(f"STARmap archive mismatch: {observed}; expected {expected}.")


def load_observed() -> tuple[pd.DataFrame, dict[str, int]]:
    with zipfile.ZipFile(RAW_ARCHIVE) as archive:
        names = set(archive.namelist())
        required = {H5AD_MEMBER, ANNOTATION_MEMBER}
        if not required.issubset(names):
            raise ValueError(f"STARmap archive lacks: {sorted(required - names)}")
        with archive.open(ANNOTATION_MEMBER) as handle:
            annotation = pd.read_csv(handle, sep="\t", index_col=0)
        with tempfile.TemporaryDirectory(prefix="simspace-starmap-") as temporary:
            extracted = Path(archive.extract(H5AD_MEMBER, path=temporary))
            with h5py.File(extracted, "r") as h5ad:
                cell_ids = decode(h5ad["obs/_index"][...])
                genes = decode(h5ad["var/_index"][...])
                layer_categories = decode(h5ad["obs/__categories/label"][...])
                layer_codes = h5ad["obs/label"][...].astype(int)
                missing_markers = sorted(set(MARKERS) - set(genes))
                if missing_markers:
                    raise ValueError(f"STARmap lacks markers: {missing_markers}")
                observed = pd.DataFrame(
                    {
                        "spatial_domain_raw": layer_categories[layer_codes],
                        "h5ad_x": h5ad["obs/X"][...],
                        "h5ad_y": h5ad["obs/Y"][...],
                    },
                    index=cell_ids,
                )
                for gene in MARKERS:
                    gene_index = int(np.flatnonzero(genes == gene)[0])
                    observed[gene] = h5ad["X"][:, gene_index]
                matrix_shape = list(h5ad["X"].shape)

    if not observed.index.equals(annotation.index):
        raise ValueError("STARmap h5ad and annotation cell orders differ.")
    observed = observed.join(annotation[["Total_counts", "X", "Y", "Annotation"]])
    if not np.allclose(observed["h5ad_x"], observed["X"]):
        raise ValueError("STARmap x coordinates disagree across files.")
    if not np.allclose(observed["h5ad_y"], observed["Y"]):
        raise ValueError("STARmap y coordinates disagree across files.")

    selected = observed[
        observed["Annotation"].isin(OBSERVED_CELL_CLASS)
        & observed["spatial_domain_raw"].isin(OBSERVED_DOMAIN)
    ].copy()
    selected["cell_class"] = selected["Annotation"].map(OBSERVED_CELL_CLASS)
    selected["spatial_domain"] = selected["spatial_domain_raw"].map(OBSERVED_DOMAIN)
    selected["x"] = selected["X"].astype(float)
    selected["y"] = -selected["Y"].astype(float)
    selected["dataset"] = "STARmap"
    counts = {
        "all_cells": len(observed),
        "selected_cells": len(selected),
        "n_genes": matrix_shape[1],
    }
    return selected, counts


def load_simulated() -> pd.DataFrame:
    cells = pd.read_csv(SIMULATED_META, sep="\t")
    expression = pd.read_csv(SIMULATED_EXPRESSION, sep="\t", usecols=list(MARKERS))
    if len(cells) != len(expression):
        raise ValueError("Figure 3 metadata and expression rows are misaligned.")
    cells = pd.concat([cells.reset_index(drop=True), expression], axis=1)
    cells["cell_class"] = cells["cell_type"].map(SIMULATED_CELL_CLASS)
    cells["spatial_domain"] = cells["niche_name"].map(SIMULATED_DOMAIN)
    if cells[["cell_class", "spatial_domain"]].isna().any().any():
        raise ValueError("An unmapped Figure 3 cell type or niche was found.")
    cells["x"] = cells["col"].astype(float)
    cells["y"] = -cells["row"].astype(float)
    cells["dataset"] = "SimSpace"
    return cells


def laminar_table(cells: pd.DataFrame, dataset: str) -> pd.DataFrame:
    counts = pd.crosstab(cells["cell_class"], cells["spatial_domain"]).reindex(
        index=CLASS_ORDER,
        columns=CLASS_ORDER,
        fill_value=0,
    )
    fractions = counts.div(counts.sum(axis=1), axis=0)
    rows = []
    for cell_class in CLASS_ORDER:
        for domain in CLASS_ORDER:
            rows.append(
                {
                    "dataset": dataset,
                    "cell_class": cell_class,
                    "spatial_domain": domain,
                    "n_cells": int(counts.loc[cell_class, domain]),
                    "within_class_fraction": float(fractions.loc[cell_class, domain]),
                }
            )
    return pd.DataFrame(rows)


def expression_table(cells: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rows = []
    for gene in MARKERS:
        gene_rows = []
        for domain in CLASS_ORDER:
            values = cells.loc[cells["spatial_domain"] == domain, gene].to_numpy(float)
            gene_rows.append(
                {
                    "dataset": dataset,
                    "gene": gene,
                    "expected_class": MARKER_EXPECTED_CLASS[gene],
                    "spatial_domain": domain,
                    "n_cells": len(values),
                    "mean_count": float(values.mean()),
                    "median_count": float(np.median(values)),
                    "detection_fraction": float(np.mean(values > 0)),
                }
            )
        denominator = sum(row["mean_count"] for row in gene_rows)
        for row in gene_rows:
            row["normalized_mean_profile"] = row["mean_count"] / denominator
        rows.extend(gene_rows)
    return pd.DataFrame(rows)


def laminar_matrix(table: pd.DataFrame, dataset: str) -> np.ndarray:
    return (
        table[table["dataset"] == dataset]
        .pivot(
            index="cell_class",
            columns="spatial_domain",
            values="within_class_fraction",
        )
        .reindex(index=CLASS_ORDER, columns=CLASS_ORDER)
        .to_numpy(float)
    )


def expression_matrix(table: pd.DataFrame, dataset: str) -> np.ndarray:
    return (
        table[table["dataset"] == dataset]
        .pivot(
            index="gene",
            columns="spatial_domain",
            values="normalized_mean_profile",
        )
        .reindex(index=MARKERS, columns=CLASS_ORDER)
        .to_numpy(float)
    )


def calculate_metrics(
    laminar: pd.DataFrame,
    expression: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    real_spatial = laminar_matrix(laminar, "STARmap")
    simulated_spatial = laminar_matrix(laminar, "SimSpace")
    real_expression = expression_matrix(expression, "STARmap")
    simulated_expression = expression_matrix(expression, "SimSpace")

    gene_rows = []
    for gene_index, gene in enumerate(MARKERS):
        observed = real_expression[gene_index]
        simulated = simulated_expression[gene_index]
        observed_peak = CLASS_ORDER[int(np.argmax(observed))]
        simulated_peak = CLASS_ORDER[int(np.argmax(simulated))]
        gene_rows.append(
            {
                "gene": gene,
                "expected_class": MARKER_EXPECTED_CLASS[gene],
                "pearson_r": float(pearsonr(observed, simulated).statistic),
                "spearman_rho": float(spearmanr(observed, simulated).statistic),
                "rmse": float(np.sqrt(np.mean((observed - simulated) ** 2))),
                "observed_peak": observed_peak,
                "simulated_peak": simulated_peak,
                "peaks_match": observed_peak == simulated_peak,
                "observed_matches_expected": observed_peak == MARKER_EXPECTED_CLASS[gene],
                "simulated_matches_expected": simulated_peak == MARKER_EXPECTED_CLASS[gene],
            }
        )
    gene_metrics = pd.DataFrame(gene_rows)

    rows = [
        {
            "comparison": "cell_class_by_spatial_domain",
            "metric": "pearson_r",
            "value": float(pearsonr(real_spatial.ravel(), simulated_spatial.ravel()).statistic),
            "n_components": 16,
        },
        {
            "comparison": "cell_class_by_spatial_domain",
            "metric": "rmse",
            "value": float(np.sqrt(np.mean((real_spatial - simulated_spatial) ** 2))),
            "n_components": 16,
        },
        {
            "comparison": "cell_class_by_spatial_domain",
            "metric": "observed_mean_diagonal_fraction",
            "value": float(np.diag(real_spatial).mean()),
            "n_components": 4,
        },
        {
            "comparison": "cell_class_by_spatial_domain",
            "metric": "simulated_mean_diagonal_fraction",
            "value": float(np.diag(simulated_spatial).mean()),
            "n_components": 4,
        },
        {
            "comparison": "eight_marker_domain_profiles",
            "metric": "pearson_r",
            "value": float(
                pearsonr(real_expression.ravel(), simulated_expression.ravel()).statistic
            ),
            "n_components": 32,
        },
        {
            "comparison": "eight_marker_domain_profiles",
            "metric": "spearman_rho",
            "value": float(
                spearmanr(real_expression.ravel(), simulated_expression.ravel()).statistic
            ),
            "n_components": 32,
        },
        {
            "comparison": "eight_marker_domain_profiles",
            "metric": "rmse",
            "value": float(np.sqrt(np.mean((real_expression - simulated_expression) ** 2))),
            "n_components": 32,
        },
        {
            "comparison": "eight_marker_domain_profiles",
            "metric": "median_gene_pearson_r",
            "value": float(gene_metrics["pearson_r"].median()),
            "n_components": len(MARKERS),
        },
        {
            "comparison": "eight_marker_domain_profiles",
            "metric": "peak_domain_match_fraction",
            "value": float(gene_metrics["peaks_match"].mean()),
            "n_components": len(MARKERS),
        },
    ]
    return pd.DataFrame(rows), gene_metrics


def scatter_classes(axis: plt.Axes, cells: pd.DataFrame, title: str) -> None:
    for cell_class in CLASS_ORDER:
        subset = cells[cells["cell_class"] == cell_class]
        axis.scatter(
            subset["x"],
            subset["y"],
            s=4.0,
            color=CLASS_COLORS[cell_class],
            edgecolors="none",
            alpha=0.88,
            rasterized=True,
        )
    axis.set_title(title, loc="left", fontweight="bold", fontsize=10.5)
    axis.set_aspect("equal")
    axis.axis("off")


def scatter_expression(
    axis: plt.Axes,
    cells: pd.DataFrame,
    gene: str,
    title: str | None = None,
):
    values = np.log1p(cells[gene].to_numpy(float))
    scale = float(np.quantile(values, 0.99))
    scaled = np.clip(values / scale, 0, 1) if scale > 0 else values
    order = np.argsort(scaled, kind="stable")
    points = axis.scatter(
        cells["x"].to_numpy()[order],
        cells["y"].to_numpy()[order],
        c=scaled[order],
        s=4.0,
        cmap=EXPRESSION_CMAP,
        vmin=0,
        vmax=1,
        edgecolors="none",
        rasterized=True,
    )
    if title is not None:
        axis.set_title(title, fontsize=9.2, pad=2, fontstyle="italic")
    axis.set_aspect("equal")
    axis.axis("off")
    return points


def annotated_heatmap(
    axis: plt.Axes,
    matrix: np.ndarray,
    row_labels: list[str],
    title: str,
    vmax: float,
    show_ylabels: bool,
) -> None:
    axis.imshow(matrix, cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            axis.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=6.6,
                color="white" if value > 0.55 * vmax else "#202020",
            )
    labels = [CLASS_LABELS[value] for value in CLASS_ORDER]
    axis.set_xticks(range(4), labels, rotation=38, ha="right", fontsize=7)
    axis.set_yticks(
        range(len(row_labels)),
        row_labels if show_ylabels else [""] * len(row_labels),
        fontsize=7,
    )
    axis.tick_params(length=0)
    axis.set_title(title, fontsize=9, pad=4)
    axis.set_xlabel("Spatial domain", fontsize=8)


def plot_gene_correlations(axis: plt.Axes, gene_metrics: pd.DataFrame) -> None:
    plot = gene_metrics.set_index("gene").reindex(MARKERS).iloc[::-1]
    colors = ["#147D92" if value >= 0 else "#D55E00" for value in plot["pearson_r"]]
    bars = axis.barh(np.arange(len(plot)), plot["pearson_r"], color=colors, height=0.65)
    axis.axvline(0, color="#303030", linewidth=0.8)
    axis.set_yticks(np.arange(len(plot)), plot.index, fontstyle="italic")
    axis.set_xlim(-1, 1)
    axis.set_xlabel("Pearson r across four domains")
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.55)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    for bar, value, matches in zip(
        bars,
        plot["pearson_r"],
        plot["peaks_match"],
        strict=True,
    ):
        offset = 0.03 if value >= 0 else -0.03
        axis.text(
            value + offset,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}" + ("  ●" if matches else ""),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=7,
        )
    axis.text(
        0.98,
        0.02,
        "● observed and simulated peak domains match",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )


def plot_figure(
    observed: pd.DataFrame,
    simulated: pd.DataFrame,
    laminar: pd.DataFrame,
    expression: pd.DataFrame,
    metrics: pd.DataFrame,
    gene_metrics: pd.DataFrame,
) -> plt.Figure:
    figure = plt.figure(figsize=(13.6, 10.8), dpi=300)
    outer_grid = figure.add_gridspec(
        3,
        1,
        height_ratios=(1.0, 0.82, 1.15),
        hspace=0.34,
    )
    top_grid = outer_grid[0].subgridspec(
        2,
        2,
        width_ratios=(1.0, 1.05),
        height_ratios=(0.12, 1.0),
        hspace=0.04,
        wspace=0.18,
    )
    spatial_panel_header = figure.add_subplot(top_grid[0, 0])
    matrix_panel_header = figure.add_subplot(top_grid[0, 1])
    for axis in (spatial_panel_header, matrix_panel_header):
        axis.axis("off")
    spatial_panel_header.text(
        0.0,
        0.5,
        "A   Spatial organization",
        transform=spatial_panel_header.transAxes,
        ha="left",
        va="center",
        fontsize=10.5,
        fontweight="bold",
    )

    spatial_map_grid = top_grid[1, 0].subgridspec(
        3,
        2,
        height_ratios=(0.09, 0.10, 1.0),
        hspace=0.01,
        wspace=0.08,
    )

    spatial_header_axis = figure.add_subplot(spatial_map_grid[0, :])
    spatial_header_axis.axis("off")
    spatial_header_axis.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=CLASS_COLORS[value],
                markeredgecolor="none",
                label=CLASS_LABELS[value],
                markersize=5,
            )
            for value in CLASS_ORDER
        ],
        loc="center left",
        bbox_to_anchor=(0.00, 0.45),
        ncol=4,
        columnspacing=0.85,
        handletextpad=0.30,
        frameon=False,
        fontsize=7.3,
        borderaxespad=0,
    )

    observed_title_axis = figure.add_subplot(spatial_map_grid[1, 0])
    simulated_title_axis = figure.add_subplot(spatial_map_grid[1, 1])
    for axis, title in (
        (observed_title_axis, "STARmap"),
        (simulated_title_axis, "SimSpace"),
    ):
        axis.axis("off")
        axis.text(
            0.0,
            0.5,
            title,
            transform=axis.transAxes,
            ha="left",
            va="center",
            fontsize=10.5,
            fontweight="bold",
        )

    observed_axis = figure.add_subplot(spatial_map_grid[2, 0])
    simulated_axis = figure.add_subplot(spatial_map_grid[2, 1])
    scatter_classes(observed_axis, observed, "")
    scatter_classes(simulated_axis, simulated, "")

    real_spatial = laminar_matrix(laminar, "STARmap")
    simulated_spatial = laminar_matrix(laminar, "SimSpace")
    spatial_r = metrics.loc[
        (metrics["comparison"] == "cell_class_by_spatial_domain")
        & (metrics["metric"] == "pearson_r"),
        "value",
    ].iloc[0]
    spatial_matrix_grid = top_grid[1, 1].subgridspec(
        2,
        2,
        height_ratios=(0.10, 1.0),
        hspace=0.02,
        wspace=0.14,
    )
    matrix_spacer_axis = figure.add_subplot(spatial_matrix_grid[0, :])
    matrix_spacer_axis.axis("off")
    real_spatial_axis = figure.add_subplot(spatial_matrix_grid[1, 0])
    simulated_spatial_axis = figure.add_subplot(spatial_matrix_grid[1, 1])
    annotated_heatmap(
        real_spatial_axis,
        real_spatial,
        [CLASS_LABELS[value] for value in CLASS_ORDER],
        "STARmap",
        1.0,
        True,
    )
    annotated_heatmap(
        simulated_spatial_axis,
        simulated_spatial,
        [CLASS_LABELS[value] for value in CLASS_ORDER],
        "SimSpace",
        1.0,
        False,
    )
    matrix_panel_header.text(
        0.0,
        0.5,
        rf"B   Cell-class distributions (matrix $r$ = {spatial_r:.3f})",
        transform=matrix_panel_header.transAxes,
        ha="left",
        va="center",
        fontweight="bold",
        fontsize=10.5,
    )

    expression_grid = outer_grid[1].subgridspec(
        1,
        2 * len(REPRESENTATIVE_MARKERS) + 1,
        width_ratios=(1.0, 1.0, 1.0, 1.0, 0.035),
        wspace=0.10,
    )
    expression_axes = []
    expression_points = None
    for column, gene in enumerate(REPRESENTATIVE_MARKERS):
        observed_expression_axis = figure.add_subplot(expression_grid[0, 2 * column])
        simulated_expression_axis = figure.add_subplot(expression_grid[0, 2 * column + 1])
        expression_points = scatter_expression(
            observed_expression_axis,
            observed,
            gene,
        )
        scatter_expression(simulated_expression_axis, simulated, gene)
        observed_expression_axis.set_title(
            rf"$\it{{{gene}}}$" + "\nSTARmap",
            fontsize=9.2,
            pad=3,
        )
        simulated_expression_axis.set_title(
            rf"$\it{{{gene}}}$" + "\nSimSpace",
            fontsize=9.2,
            pad=3,
        )
        expression_axes.extend([observed_expression_axis, simulated_expression_axis])

    expression_axes[0].text(
        0.0,
        1.23,
        "C   Representative marker expression",
        transform=expression_axes[0].transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    expression_colorbar_axis = figure.add_subplot(expression_grid[0, -1])
    colorbar = figure.colorbar(
        expression_points,
        cax=expression_colorbar_axis,
        orientation="vertical",
    )
    colorbar.set_label(
        "Scaled log1p count\n(within dataset/gene)",
        fontsize=7.5,
        rotation=270,
        labelpad=19,
    )
    colorbar.ax.tick_params(labelsize=7, length=2)
    colorbar.set_ticks([0, 0.5, 1])

    real_expression = expression_matrix(expression, "STARmap")
    simulated_expression = expression_matrix(expression, "SimSpace")
    expression_r = metrics.loc[
        (metrics["comparison"] == "eight_marker_domain_profiles")
        & (metrics["metric"] == "pearson_r"),
        "value",
    ].iloc[0]
    summary_grid = outer_grid[2].subgridspec(
        2,
        2,
        width_ratios=(1.0, 1.10),
        height_ratios=(0.12, 1.0),
        hspace=0.05,
        wspace=0.24,
    )
    profile_panel_header = figure.add_subplot(summary_grid[0, 0])
    correlation_panel_header = figure.add_subplot(summary_grid[0, 1])
    for axis in (profile_panel_header, correlation_panel_header):
        axis.axis("off")
    profile_panel_header.text(
        0.0,
        0.5,
        rf"D   All eight shared marker profiles (matrix $r$ = {expression_r:.3f})",
        transform=profile_panel_header.transAxes,
        ha="left",
        va="center",
        fontweight="bold",
        fontsize=10.5,
    )
    correlation_panel_header.text(
        0.0,
        0.5,
        "E   Per-marker profile agreement",
        transform=correlation_panel_header.transAxes,
        ha="left",
        va="center",
        fontweight="bold",
        fontsize=10.5,
    )

    profile_grid = summary_grid[1, 0].subgridspec(1, 2, wspace=0.12)
    real_expression_axis = figure.add_subplot(profile_grid[0, 0])
    simulated_expression_axis = figure.add_subplot(profile_grid[0, 1])
    annotated_heatmap(
        real_expression_axis,
        real_expression,
        list(MARKERS),
        "STARmap",
        0.65,
        True,
    )
    annotated_heatmap(
        simulated_expression_axis,
        simulated_expression,
        list(MARKERS),
        "SimSpace",
        0.65,
        False,
    )
    correlation_axis = figure.add_subplot(summary_grid[1, 1])
    plot_gene_correlations(correlation_axis, gene_metrics)

    figure.subplots_adjust(left=0.050, right=0.965, bottom=0.060, top=0.965)
    summary_shift = 0.025
    for axis in (
        profile_panel_header,
        correlation_panel_header,
        real_expression_axis,
        simulated_expression_axis,
        correlation_axis,
    ):
        position = axis.get_position()
        axis.set_position(
            [
                position.x0,
                position.y0 + summary_shift,
                position.width,
                position.height,
            ]
        )
    return figure


def write_outputs(
    observed: pd.DataFrame,
    simulated: pd.DataFrame,
    observed_counts: dict[str, int],
    laminar: pd.DataFrame,
    expression: pd.DataFrame,
    metrics: pd.DataFrame,
    gene_metrics: pd.DataFrame,
) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    observed_columns = [
        "Annotation",
        "spatial_domain_raw",
        "cell_class",
        "spatial_domain",
        "x",
        "y",
        *MARKERS,
    ]
    observed[observed_columns].to_csv(
        DATA_DIR / "observed_mapped_cells.tsv.gz",
        sep="\t",
        index=True,
        index_label="cell_id",
        compression="gzip",
    )
    simulated_columns = [
        "cell_type",
        "niche_name",
        "cell_class",
        "spatial_domain",
        "x",
        "y",
        *MARKERS,
    ]
    simulated[simulated_columns].to_csv(
        DATA_DIR / "simulated_mapped_cells.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    laminar.to_csv(DATA_DIR / "laminar_enrichment.tsv", sep="\t", index=False)
    expression.to_csv(DATA_DIR / "marker_profiles.tsv", sep="\t", index=False)
    metrics.to_csv(DATA_DIR / "comparison_metrics.tsv", sep="\t", index=False)
    gene_metrics.to_csv(DATA_DIR / "gene_profile_metrics.tsv", sep="\t", index=False)

    configuration = {
        "status": "selected Supplementary Figure 6 analysis",
        "reviewer_comment": "R2-M4",
        "analysis_role": (
            "Post-generation comparison of the fixed reference-free Figure 3 "
            "anchor; STARmap was not used to fit, tune, or select the simulated anchor."
        ),
        "observed_dataset": {
            "name": "STARmap mouse visual cortex BY3 1k-gene sample",
            "zenodo_record": SOURCE_RECORD,
            "source_url": SOURCE_URL,
            "zenodo_doi": SOURCE_DOI,
            "original_starmap_doi": ORIGINAL_STARMAP_DOI,
            "license": SOURCE_LICENSE,
            "archive_bytes": EXPECTED_RAW_BYTES,
            "archive_sha256": EXPECTED_RAW_SHA256,
            "archive_md5": EXPECTED_RAW_MD5,
            "matrix_cells_by_genes": [observed_counts["all_cells"], observed_counts["n_genes"]],
            "selected_cells": observed_counts["selected_cells"],
        },
        "selection_rule": (
            "Cells with an excitatory layer annotation eL2/3, eL4, eL5, "
            "eL6-1, or eL6-2 and a cortical domain label L2/3, L4, L5, or L6."
        ),
        "class_order": list(CLASS_ORDER),
        "observed_cell_class_mapping": OBSERVED_CELL_CLASS,
        "observed_spatial_domain_mapping": OBSERVED_DOMAIN,
        "simulated_cell_class_mapping": SIMULATED_CELL_CLASS,
        "simulated_spatial_domain_mapping": SIMULATED_DOMAIN,
        "markers": list(MARKERS),
        "representative_map_markers": list(REPRESENTATIVE_MARKERS),
        "representative_map_note": (
            "Cux2 and Rorb are displayed as spatial examples; all eight shared "
            "markers remain in the profile panels and numeric metrics."
        ),
        "marker_expected_class": MARKER_EXPECTED_CLASS,
        "expression_comparison": (
            "For each gene, domain-wise mean counts are normalized within "
            "each dataset to sum to one. Maps show log1p counts scaled to "
            "each dataset and gene's 99th percentile."
        ),
        "independence_note": (
            "One observed sample and one fixed simulated realization; cells "
            "and domains are not treated as biological replicates."
        ),
    }
    (DATA_DIR / "analysis_config.json").write_text(
        json.dumps(configuration, indent=2) + "\n",
        encoding="utf-8",
    )

    manifest = []
    for role, path, source in (
        ("observed_archive", RAW_ARCHIVE, SOURCE_URL),
        ("simulated_metadata", SIMULATED_META, "main_figures/Fig3"),
        ("simulated_expression", SIMULATED_EXPRESSION, "main_figures/Fig3"),
    ):
        manifest.append(
            {
                "role": role,
                "path": str(path.relative_to(REPO_ROOT)),
                "source": source,
                "bytes": path.stat().st_size,
                "sha256": checksum(path),
            }
        )
    pd.DataFrame(manifest).to_csv(
        DATA_DIR / "input_manifest.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {"software": "Python", "version": platform.python_version()},
            {"software": "h5py", "version": h5py.__version__},
            {"software": "NumPy", "version": np.__version__},
            {"software": "pandas", "version": pd.__version__},
            {"software": "SciPy", "version": scipy.__version__},
            {"software": "matplotlib", "version": matplotlib.__version__},
            {"software": "platform", "version": platform.platform()},
            {"software": "python_executable", "version": sys.executable},
        ]
    ).to_csv(DATA_DIR / "software_versions.tsv", sep="\t", index=False)


def main() -> None:
    verify_inputs()
    observed, observed_counts = load_observed()
    simulated = load_simulated()
    laminar = pd.concat(
        [laminar_table(observed, "STARmap"), laminar_table(simulated, "SimSpace")],
        ignore_index=True,
    )
    expression = pd.concat(
        [expression_table(observed, "STARmap"), expression_table(simulated, "SimSpace")],
        ignore_index=True,
    )
    metrics, gene_metrics = calculate_metrics(laminar, expression)
    write_outputs(
        observed,
        simulated,
        observed_counts,
        laminar,
        expression,
        metrics,
        gene_metrics,
    )

    figure = plot_figure(observed, simulated, laminar, expression, metrics, gene_metrics)
    output = SCRIPT_DIR / "SFig6.png"
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(output, EXAMPLE_DIR / output.name)

    lookup = metrics.set_index(["comparison", "metric"])["value"]
    print(
        f"Generated {output.name} from {len(observed):,} STARmap and "
        f"{len(simulated):,} simulated cells; spatial matrix r="
        f"{lookup[('cell_class_by_spatial_domain', 'pearson_r')]:.3f}; "
        f"eight-marker matrix r="
        f"{lookup[('eight_marker_domain_profiles', 'pearson_r')]:.3f}; "
        f"peak matches={gene_metrics['peaks_match'].sum()}/{len(gene_metrics)}."
    )


if __name__ == "__main__":
    main()
