#!/usr/bin/env python3
"""Generate the three standalone panels for Figure 3.

Panel A shows the four reference-free niche bands used to construct the
layered cortex. Panel B varies one geometry parameter and one biological
organization parameter in a 5 x 5 design whose center is the selected anchor.
Panel C attaches a reference-free, cortex-inspired molecular profile to the
anchor and displays two canonical layer-enriched marker genes.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from cortex_reference_free import (
    CELL_COLORS,
    CELL_TYPES,
    DEFAULT_NICHE_MIXTURE,
    CortexConfig,
    plot_niches,
    ss,
    simulate_cortex,
    summarize_simulation,
)


INNER_RADIUS_VALUES = np.asarray([0.58, 0.51, 0.44, 0.37, 0.30])
SPECIFICITY_VALUES = np.asarray([3.00, 1.50, 1.00, 0.50, 0.00])
ANCHOR_GEOMETRY_INDEX = 2
ANCHOR_SPECIFICITY_INDEX = 2

MARKER_TARGETS = (
    ("Cux2", "L2/3 IT"),
    ("Rorb", "L4/5 IT"),
    ("Fezf2", "L5 ET"),
    ("Deptor", "L5 IT"),
    ("Tshz2", "L5/6 NP"),
    ("Foxp2", "L6 CT"),
    ("Sulf2", "L6 IT"),
    ("Car3", "L6 IT Car3"),
    ("Ctgf", "L6b"),
)

# Relative expression patterns are deliberately qualitative. They encode
# well-established subclass enrichment and limited expression in related
# cortical subclasses; they are not estimates fitted from a count matrix.
MARKER_WEIGHTS = {
    "Cux2": [1.00, 0.45, 0.04, 0.04, 0.02, 0.01, 0.01, 0.01, 0.01],
    "Rorb": [0.20, 1.00, 0.25, 0.45, 0.08, 0.02, 0.03, 0.01, 0.01],
    "Fezf2": [0.01, 0.04, 1.00, 0.15, 0.16, 0.12, 0.04, 0.02, 0.02],
    "Deptor": [0.01, 0.06, 0.20, 1.00, 0.18, 0.08, 0.16, 0.04, 0.02],
    "Tshz2": [0.01, 0.02, 0.10, 0.08, 1.00, 0.20, 0.10, 0.04, 0.08],
    "Foxp2": [0.005, 0.01, 0.06, 0.04, 0.16, 1.00, 0.15, 0.03, 0.25],
    "Sulf2": [0.005, 0.01, 0.04, 0.12, 0.12, 0.15, 1.00, 0.18, 0.12],
    "Car3": [0.005, 0.005, 0.02, 0.03, 0.08, 0.04, 0.12, 1.00, 0.08],
    "Ctgf": [0.005, 0.005, 0.01, 0.02, 0.08, 0.15, 0.08, 0.06, 1.00],
}

EXPRESSION_CMAP = LinearSegmentedColormap.from_list(
    "cortex_expression",
    ["#e6e8eb", "#d5e7f1", "#9ecae1", "#4292c6", "#08519c"],
)


def _set_cortex_limits(ax: plt.Axes, config: CortexConfig) -> None:
    ax.set_xlim(-4, config.cols + 3)
    ax.set_ylim(-config.rows - 3, 4)
    ax.set_aspect("equal")
    ax.axis("off")


def _scatter_cell_types(
    ax: plt.Axes,
    simulation,
    config: CortexConfig,
    point_size: float,
) -> None:
    state = simulation.meta["state"].astype(int).to_numpy()
    x = simulation.meta["col"].to_numpy()
    y = -simulation.meta["row"].to_numpy()
    for index, color in enumerate(CELL_COLORS):
        selected = state == index
        ax.scatter(
            x[selected],
            y[selected],
            s=point_size,
            color=color,
            edgecolors="none",
            alpha=0.94,
            rasterized=True,
        )
    _set_cortex_limits(ax, config)


def layer_specificity_mixture(factor: float) -> np.ndarray:
    """Interpolate/extrapolate niche mixtures around the selected anchor."""

    if np.isclose(factor, 1.0):
        return DEFAULT_NICHE_MIXTURE.copy()
    global_profile = DEFAULT_NICHE_MIXTURE.mean(axis=0, keepdims=True)
    mixture = global_profile + factor * (DEFAULT_NICHE_MIXTURE - global_profile)
    mixture = np.clip(mixture, 1e-5, None)
    return mixture / mixture.sum(axis=1, keepdims=True)


def make_panel_b(
    base_config: CortexConfig,
    anchor_simulation,
    output_path: Path,
    data_dir: Path,
    dpi: int,
) -> pd.DataFrame:
    """Draw a 5 x 5 geometry-by-layer-specificity parameter sweep."""

    fig, axes = plt.subplots(5, 5, figsize=(11.4, 10.9), dpi=190)
    records: list[dict[str, float | int | bool]] = []

    for row_index, specificity in enumerate(SPECIFICITY_VALUES):
        mixture = layer_specificity_mixture(float(specificity))
        for col_index, inner_radius in enumerate(INNER_RADIUS_VALUES):
            config = replace(base_config, inner_radius=float(inner_radius))
            is_anchor = (
                row_index == ANCHOR_SPECIFICITY_INDEX
                and col_index == ANCHOR_GEOMETRY_INDEX
            )
            if is_anchor:
                simulation = anchor_simulation
            else:
                simulation, _, _ = simulate_cortex(config, mixture)
            _, diagnostics = summarize_simulation(simulation)
            ax = axes[row_index, col_index]
            _scatter_cell_types(ax, simulation, config, point_size=1.75)

            if is_anchor:
                ax.add_patch(
                    Rectangle(
                        (0.006, 0.006),
                        0.988,
                        0.988,
                        transform=ax.transAxes,
                        fill=False,
                        edgecolor="#202124",
                        linewidth=1.5,
                        clip_on=False,
                        zorder=20,
                    )
                )
                ax.text(
                    0.04,
                    0.95,
                    "anchor",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=10.0,
                    fontweight="bold",
                    color="#202124",
                    bbox={
                        "boxstyle": "round,pad=0.20",
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.90,
                    },
                    zorder=21,
                )

            if row_index == 0:
                ax.set_title(f"{inner_radius:.2f}", fontsize=10, pad=3.0)
            if col_index == 0:
                ax.text(
                    -0.08,
                    0.50,
                    f"{specificity:.1f}×",
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=10,
                )

            records.append(
                {
                    "row": row_index + 1,
                    "column": col_index + 1,
                    "inner_radius": float(inner_radius),
                    "layer_specificity_factor": float(specificity),
                    "is_anchor": is_anchor,
                    "n_cells": diagnostics["n_cells"],
                    "composition_rmse_visual_guide": diagnostics[
                        "composition_rmse_visual_guide"
                    ],
                    "median_depth_mae_visual_guide": diagnostics[
                        "median_depth_mae_visual_guide"
                    ],
                    "laminar_ordering_violations": diagnostics[
                        "laminar_ordering_violations"
                    ],
                }
            )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="none",
            markersize=7,
            label=cell_type,
        )
        for cell_type, color in zip(CELL_TYPES, CELL_COLORS)
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.015),
        ncol=5,
        frameon=False,
        fontsize=10,
        handletextpad=0.35,
        columnspacing=1.25,
    )
    fig.supxlabel(
        "Cortical thickness (thin → thick; inner-radius parameter)",
        x=0.535,
        y=0.083,
        fontsize=14,
    )
    fig.supylabel(
        "Layer specificity (low → high)",
        x=0.018,
        y=0.535,
        fontsize=14,
    )
    fig.subplots_adjust(
        left=0.085,
        right=0.995,
        top=0.965,
        bottom=0.125,
        wspace=0.015,
        hspace=0.045,
    )
    fig.savefig(output_path, dpi=dpi, facecolor="white")
    plt.close(fig)

    sweep = pd.DataFrame.from_records(records)
    data_dir.mkdir(parents=True, exist_ok=True)
    sweep.to_csv(data_dir / "sweep_parameters.tsv", sep="\t", index=False)
    return sweep


def build_gene_parameters(seed: int = 7301) -> pd.DataFrame:
    """Create a cortex-inspired Gamma-Poisson mean-expression design."""

    rng = np.random.default_rng(seed)
    cell_index = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    genes: list[tuple[str, int, np.ndarray | None, str]] = []

    for gene, target in MARKER_TARGETS:
        genes.append(
            (
                gene,
                cell_index[target],
                np.asarray(MARKER_WEIGHTS[gene], dtype=float),
                "named cortical subclass marker",
            )
        )

    for target_index, cell_type in enumerate(CELL_TYPES):
        for replicate in range(1, 3):
            weights = np.full(len(CELL_TYPES), 0.025, dtype=float)
            weights[target_index] = 1.0
            genes.append(
                (
                    f"Marker_{target_index + 1:02d}_{replicate}",
                    target_index,
                    weights,
                    f"generic marker for {cell_type}",
                )
            )

    for background_index in range(1, 73):
        genes.append(
            (
                f"Background_{background_index:03d}",
                -1,
                None,
                "background gene",
            )
        )

    rows: list[dict[str, float | int | str]] = []
    for gene_id, (gene, target, weights, role) in enumerate(genes):
        baseline = rng.gamma(shape=1.45, scale=0.48, size=len(CELL_TYPES)) + 0.04
        if weights is not None:
            marker_amplitude = float(rng.gamma(shape=8.0, scale=1.0))
            means = baseline + marker_amplitude * weights
        else:
            marker_amplitude = 0.0
            means = baseline
        row: dict[str, float | int | str] = {
            "GeneID": gene_id,
            "gene": gene,
            "Marker": target,
            "LRindex": -1,
            "role": role,
            "marker_amplitude": marker_amplitude,
        }
        for state_index, mean in enumerate(means):
            row[f"Type_{state_index}"] = float(mean)
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def simulate_molecular_profile(
    simulation,
    output_dir: Path,
    parameter_seed: int = 7301,
    count_seed: int = 7302,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach reference-free expression to the selected spatial realization."""

    gene_parameters = build_gene_parameters(seed=parameter_seed)
    simspace_columns = [
        "GeneID",
        "Marker",
        "LRindex",
        *[f"Type_{index}" for index in range(len(CELL_TYPES))],
    ]
    simulation.gene_meta = gene_parameters.copy()

    # SimSpace's reference-free molecular engine draws Poisson counts from the
    # cell-type-conditional mean matrix. The means above were themselves drawn
    # from Gamma distributions, yielding a declared Gamma-Poisson design.
    simulation_omics = ss.omics.simOmics(
        gene_parameters[simspace_columns],
        simulation.meta,
        seed=count_seed,
    )
    rename = {
        f"Gene_{gene_id}": gene
        for gene_id, gene in zip(gene_parameters["GeneID"], gene_parameters["gene"])
    }
    counts = simulation_omics.rename(columns=rename).astype(np.int32)
    simulation.omics = counts

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_parameters.to_csv(
        output_dir / "gene_parameters.tsv",
        sep="\t",
        index=False,
    )
    counts.to_csv(
        output_dir / "simulated_expression.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    marker_names = [gene for gene, _ in MARKER_TARGETS]
    marker_summary = (
        counts[marker_names]
        .assign(cell_type=simulation.meta["cell_type"].astype(str).to_numpy())
        .groupby("cell_type", sort=False)[marker_names]
        .mean()
        .reindex(CELL_TYPES)
    )
    marker_summary.index.name = "cell_type"
    marker_summary.to_csv(output_dir / "marker_expression_by_cell_type.tsv", sep="\t")
    return counts, gene_parameters


def make_panel_c(
    simulation,
    counts: pd.DataFrame,
    output_path: Path,
    config: CortexConfig,
    dpi: int,
) -> None:
    """Plot superficial and deep marker expression on the anchor cortex."""

    display_genes = (("Cux2", "L2/3 IT-enriched"), ("Foxp2", "L6 CT-enriched"))
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.35), dpi=210)
    x = simulation.meta["col"].to_numpy()
    y = -simulation.meta["row"].to_numpy()

    for ax, (gene, description) in zip(axes, display_genes):
        values = counts[gene].to_numpy(dtype=float)
        order = np.argsort(values, kind="stable")
        vmax = max(4.0, float(np.quantile(values, 0.995)))
        scatter = ax.scatter(
            x[order],
            y[order],
            c=values[order],
            s=5.6,
            cmap=EXPRESSION_CMAP,
            norm=Normalize(vmin=0.0, vmax=vmax),
            edgecolors="none",
            rasterized=True,
        )
        _set_cortex_limits(ax, config)
        ax.set_title(
            rf"$\it{{{gene}}}$" + f"\n{description} marker",
            fontsize=14,
            pad=5.0,
        )
        colorbar = fig.colorbar(scatter, ax=ax, fraction=0.043, pad=0.018)
        colorbar.set_label("Simulated count", fontsize=10)
        colorbar.ax.tick_params(labelsize=10, length=2.2)
        colorbar.outline.set_linewidth(0.5)

    fig.subplots_adjust(left=0.005, right=0.985, top=0.88, bottom=0.015, wspace=0.10)
    fig.savefig(output_path, dpi=dpi, facecolor="white", bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)


def write_figure_brief(output_path: Path, config: CortexConfig) -> None:
    payload = {
        "figure_type": "single-cell spatial simulation demonstration",
        "manuscript_role": "main Figure 3",
        "one_sentence_takeaway": (
            "Reference-free SimSpace parameters generate a designed layered cortex, "
            "a controllable family of related tissues, and spatially organized "
            "cell-type-conditional molecular profiles."
        ),
        "panels": {
            "A": "Four de-novo cortical niches in the selected anchor realization.",
            "B": (
                "A 5 x 5 sweep of cortical thickness and layer specificity; the center "
                "is the selected anchor and all other random seeds are held fixed."
            ),
            "C": (
                "Reference-free Gamma-Poisson expression of Cux2 and Foxp2 on the "
                "selected anchor, conditional on the simulated cell labels."
            ),
        },
        "spatial_reference_used_for_generation": False,
        "expression_reference_used_for_generation": False,
        "anchor_config": asdict(config),
        "inner_radius_values_thin_to_thick": INNER_RADIUS_VALUES.tolist(),
        "layer_specificity_values_top_to_bottom": SPECIFICITY_VALUES.tolist(),
        "molecular_assumption": (
            "Named markers encode qualitative cortical-subclass enrichment and "
            "are not quantitative estimates fitted from mouse scRNA-seq."
        ),
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for Panel_A.png, Panel_B.png, Panel_C.png, and panel data.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config = CortexConfig()
    anchor, _, _ = simulate_cortex(config)

    plot_niches(anchor, args.output_dir / "Panel_A.png")
    anchor_data_dir = args.output_dir / "Panel_A_C_data"
    anchor_data_dir.mkdir(parents=True, exist_ok=True)
    anchor.meta.to_csv(
        anchor_data_dir / "anchor_cells.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    sweep = make_panel_b(
        config,
        anchor,
        args.output_dir / "Panel_B.png",
        args.output_dir / "Panel_B_data",
        dpi=args.dpi,
    )
    counts, _ = simulate_molecular_profile(
        anchor,
        args.output_dir / "Panel_C_data",
    )
    make_panel_c(
        anchor,
        counts,
        args.output_dir / "Panel_C.png",
        config,
        dpi=args.dpi,
    )
    write_figure_brief(args.output_dir / "figure_brief.json", config)

    anchor_row = sweep.loc[sweep["is_anchor"]].iloc[0]
    print(
        "Generated Panel_A.png, Panel_B.png, and Panel_C.png; "
        f"anchor n={int(anchor_row['n_cells']):,}."
    )


if __name__ == "__main__":
    main()
