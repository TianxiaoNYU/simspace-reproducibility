#!/usr/bin/env python3
"""Reference-free SimSpace experiment for a layered mouse-cortex-like tissue.

The spatial coordinates and four curved niches are generated de novo from a
small set of geometric parameters. Nine cortical excitatory neuron classes are
then sampled and refined with SimSpace's niche-specific MRF.

This is an experiment driver, not a spatial-reference fitting workflow. The
optional target screenshot is used only when drawing a visual comparison.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter


WORKSPACE = Path(__file__).resolve().parents[3]
LOCAL_SIMSPACE = WORKSPACE / "SimSpace"
if str(LOCAL_SIMSPACE) not in sys.path:
    sys.path.insert(0, str(LOCAL_SIMSPACE))

import simspace as ss  # noqa: E402


CELL_TYPES = (
    "L2/3 IT",
    "L4/5 IT",
    "L5 ET",
    "L5 IT",
    "L5/6 NP",
    "L6 CT",
    "L6 IT",
    "L6 IT Car3",
    "L6b",
)

# Colors sampled from the supplied scHolography Figure 3A screenshot.
CELL_COLORS = (
    "#bfd6e6",
    "#4579b1",
    "#dbee9b",
    "#56a539",
    "#f3acac",
    "#c73b2a",
    "#f7cb80",
    "#e28030",
    "#dec3ee",
)

NICHE_NAMES = (
    "superficial (L2/3)",
    "middle (L4/5)",
    "deep (L5)",
    "inner (L5/6-L6b)",
)

# Rows are niches; columns follow CELL_TYPES. These are generative priors, not
# proportions measured from a spatial reference.
DEFAULT_NICHE_MIXTURE = np.asarray(
    [
        [0.782, 0.15, 0.025, 0.025, 0.008, 0.005, 0.003, 0.002, 0.000],
        [0.20, 0.505, 0.085, 0.12, 0.040, 0.020, 0.015, 0.010, 0.000],
        [0.025, 0.10, 0.25, 0.218, 0.20, 0.140, 0.040, 0.025, 0.002],
        [0.000, 0.015, 0.030, 0.060, 0.33, 0.340, 0.120, 0.065, 0.040],
    ],
    dtype=float,
)
DEFAULT_NICHE_MIXTURE /= DEFAULT_NICHE_MIXTURE.sum(axis=1, keepdims=True)

# An approximate visual target used only for tuning diagnostics. It is not
# presented as a quantitative extraction from the published panel.
VISUAL_COMPOSITION_GUIDE = np.asarray(
    [0.31, 0.18, 0.10, 0.11, 0.12, 0.09, 0.045, 0.025, 0.01],
    dtype=float,
)
VISUAL_COMPOSITION_GUIDE /= VISUAL_COMPOSITION_GUIDE.sum()

EXPECTED_DEPTH = np.asarray([0.14, 0.34, 0.50, 0.54, 0.68, 0.80, 0.82, 0.86, 0.94])
LAMINAR_RANK = np.asarray([0, 1, 2, 2, 3, 4, 4, 4, 5], dtype=float)


@dataclass(frozen=True)
class CortexConfig:
    """Parameters controlling geometry, niches, cell MRF, and sampling."""

    rows: int = 104
    cols: int = 118
    center_x: float = 1.15
    center_y: float = 1.18
    y_scale: float = 1.05
    outer_radius_at_top: float = 1.14
    outer_radius_slope: float = 0.18
    inner_radius: float = 0.44
    angle_min_deg: float = -169.0
    angle_max_deg: float = -82.0
    niche_cut_1: float = 0.29
    niche_cut_2: float = 0.51
    niche_cut_3: float = 0.69
    boundary_noise: float = 0.010
    interface_noise: float = 0.018
    geometry_seed: int = 419
    simulation_seed: int = 23
    initial_seed: int = 1023
    neighbor_distance: int = 2
    num_iterations: int = 1
    phi: float = 2.60
    theta_bias_strength: float = 0.20
    theta_depth_penalty: float = 0.12
    rare_state_cutoff: float = 0.015
    rare_state_penalty: float = 0.40
    retention: float = 0.34
    coordinate_jitter: float = 0.22


def _smooth_noise(shape: tuple[int, int], seed: int, sigma: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    field = gaussian_filter(rng.normal(size=shape), sigma=sigma, mode="reflect")
    std = float(field.std())
    return field / std if std > 0 else field


def build_cortex_geometry(
    config: CortexConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return tissue mask, four niche labels, and normalized cortical depth."""

    row, col = np.indices((config.rows, config.cols))
    x = col / (config.cols - 1)
    y = row / (config.rows - 1)
    dx = x - config.center_x
    dy = (y - config.center_y) / config.y_scale
    radius = np.sqrt(dx**2 + dy**2)
    angle = np.degrees(np.arctan2(dy, dx))

    angle_span = config.angle_max_deg - config.angle_min_deg
    angular_position = np.clip(
        (config.angle_max_deg - angle) / angle_span,
        0.0,
        1.0,
    )
    outer_radius = (
        config.outer_radius_at_top
        + config.outer_radius_slope * angular_position
        + 0.012 * np.sin(np.radians(3.0 * angle + 25.0))
    )
    boundary_field = _smooth_noise(
        (config.rows, config.cols), config.geometry_seed, sigma=5.0
    )
    outer_radius = outer_radius + config.boundary_noise * boundary_field

    tissue = (
        (angle >= config.angle_min_deg)
        & (angle <= config.angle_max_deg)
        & (radius >= config.inner_radius)
        & (radius <= outer_radius)
    )
    depth = np.clip(
        (outer_radius - radius) / (outer_radius - config.inner_radius),
        0.0,
        1.0,
    )

    interface_field = _smooth_noise(
        (config.rows, config.cols), config.geometry_seed + 1, sigma=7.0
    )
    interface_wave = (
        config.interface_noise * interface_field
        + 0.012 * np.sin(np.radians(4.0 * angle - 10.0))
    )
    cuts = np.stack(
        [
            np.full_like(depth, config.niche_cut_1) + interface_wave,
            np.full_like(depth, config.niche_cut_2) + 0.65 * interface_wave,
            np.full_like(depth, config.niche_cut_3) + 0.45 * interface_wave,
        ],
        axis=0,
    )
    cuts = np.sort(cuts, axis=0)
    niches = np.zeros_like(row, dtype=int)
    niches[depth >= cuts[0]] = 1
    niches[depth >= cuts[1]] = 2
    niches[depth >= cuts[2]] = 3
    niches[~tissue] = 0
    depth[~tissue] = np.nan
    return tissue, niches, depth


def build_theta_matrices(
    niche_mixture: np.ndarray,
    bias_strength: float,
    depth_penalty: float,
    rare_state_cutoff: float,
    rare_state_penalty: float,
) -> list[np.ndarray]:
    """Construct symmetric niche-specific phenotype affinity matrices."""

    matrices: list[np.ndarray] = []
    for probabilities in niche_mixture:
        log_probability = np.log(np.clip(probabilities, 1e-5, None))
        centered = log_probability - log_probability.mean()
        scale = np.max(np.abs(centered))
        abundance_bias = centered / scale if scale > 0 else centered
        theta = np.zeros((len(CELL_TYPES), len(CELL_TYPES)), dtype=float)
        for i in range(len(CELL_TYPES)):
            for j in range(i + 1, len(CELL_TYPES)):
                laminar_distance = abs(LAMINAR_RANK[i] - LAMINAR_RANK[j])
                value = (
                    0.10
                    - depth_penalty * laminar_distance
                    + bias_strength * (abundance_bias[i] + abundance_bias[j]) / 2.0
                )
                if (
                    probabilities[i] < rare_state_cutoff
                    or probabilities[j] < rare_state_cutoff
                ):
                    value -= rare_state_penalty
                theta[i, j] = theta[j, i] = np.clip(value, -0.80, 0.65)
        np.fill_diagonal(theta, 1.0)
        matrices.append(theta)
    return matrices


def simulate_cortex(
    config: CortexConfig,
    niche_mixture: np.ndarray = DEFAULT_NICHE_MIXTURE,
) -> tuple[ss.SimSpace, np.ndarray, np.ndarray]:
    """Generate one de-novo cortex realization using SimSpace."""

    tissue, niches, depth = build_cortex_geometry(config)
    theta = build_theta_matrices(
        niche_mixture,
        bias_strength=config.theta_bias_strength,
        depth_penalty=config.theta_depth_penalty,
        rare_state_cutoff=config.rare_state_cutoff,
        rare_state_penalty=config.rare_state_penalty,
    )
    simulation = ss.SimSpace(
        shape=(config.rows, config.cols),
        num_states=len(CELL_TYPES),
        num_iterations=config.num_iterations,
        theta=theta,
        phi=config.phi,
        neighborhood=ss.spatial.generate_offsets(
            config.neighbor_distance, method="manhattan"
        ),
        random_seed=config.simulation_seed,
    )

    # A niche-dependent categorical initialization is part of the declared
    # reference-free parameter set. SimSpace's MRF then imposes local cell-cell
    # dependence within each niche.
    rng = np.random.default_rng(config.initial_seed)
    grid = np.full((config.rows, config.cols), -1, dtype=int)
    for niche_index in range(len(NICHE_NAMES)):
        locations = np.flatnonzero(tissue & (niches == niche_index))
        grid.flat[locations] = rng.choice(
            len(CELL_TYPES),
            size=len(locations),
            p=niche_mixture[niche_index],
        )
    simulation.grid = grid
    simulation.niche = niches
    simulation.gibbs_sampler()
    simulation.density_sampler([float(config.retention)] * len(CELL_TYPES))
    simulation.perturbation(step=config.coordinate_jitter)

    lattice_row = np.clip(np.rint(simulation.meta["row"]).astype(int), 0, config.rows - 1)
    lattice_col = np.clip(np.rint(simulation.meta["col"]).astype(int), 0, config.cols - 1)
    state_int = simulation.meta["state"].astype(int)
    simulation.meta["cell_type"] = pd.Categorical(
        [CELL_TYPES[index] for index in state_int],
        categories=CELL_TYPES,
        ordered=True,
    )
    simulation.meta["niche_name"] = pd.Categorical(
        [NICHE_NAMES[index] for index in simulation.meta["niche"].astype(int)],
        categories=NICHE_NAMES,
        ordered=True,
    )
    simulation.meta["cortical_depth"] = depth[lattice_row, lattice_col]
    return simulation, tissue, depth


def summarize_simulation(simulation: ss.SimSpace) -> tuple[pd.DataFrame, dict[str, float]]:
    counts = (
        simulation.meta["cell_type"]
        .value_counts(sort=False)
        .reindex(CELL_TYPES, fill_value=0)
    )
    proportions = counts / counts.sum()
    median_depth = (
        simulation.meta.groupby("cell_type", observed=False)["cortical_depth"]
        .median()
        .reindex(CELL_TYPES)
    )
    composition_rmse = float(
        np.sqrt(np.mean((proportions.to_numpy() - VISUAL_COMPOSITION_GUIDE) ** 2))
    )
    depth_mae = float(np.nanmean(np.abs(median_depth.to_numpy() - EXPECTED_DEPTH)))
    ordering_violations = int(
        np.sum(np.diff(median_depth.to_numpy()[[0, 1, 2, 4, 5, 8]]) < 0)
    )
    diagnostics = {
        "n_cells": int(counts.sum()),
        "composition_rmse_visual_guide": composition_rmse,
        "median_depth_mae_visual_guide": depth_mae,
        "laminar_ordering_violations": ordering_violations,
        "tuning_score": composition_rmse + 0.45 * depth_mae + 0.03 * ordering_violations,
    }
    table = pd.DataFrame(
        {
            "cell_type": CELL_TYPES,
            "count": counts.to_numpy(),
            "proportion": proportions.to_numpy(),
            "visual_composition_guide": VISUAL_COMPOSITION_GUIDE,
            "median_cortical_depth": median_depth.to_numpy(),
            "expected_depth_guide": EXPECTED_DEPTH,
        }
    )
    return table, diagnostics


def _scatter_simulation(ax: plt.Axes, simulation: ss.SimSpace, point_size: float = 5.2) -> None:
    for index, (cell_type, color) in enumerate(zip(CELL_TYPES, CELL_COLORS)):
        subset = simulation.meta[simulation.meta["cell_type"] == cell_type]
        ax.scatter(
            subset["col"],
            -subset["row"],
            s=point_size,
            color=color,
            edgecolors="none",
            alpha=0.92,
            rasterized=True,
            label=cell_type,
        )
    ax.set_aspect("equal")
    ax.axis("off")


def plot_single(
    simulation: ss.SimSpace,
    output_path: Path,
    title: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.5), dpi=240)
    _scatter_simulation(ax, simulation)
    if title:
        ax.set_title(title, fontsize=11)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        markerscale=2.2,
        fontsize=8,
        handletextpad=0.35,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_niches(simulation: ss.SimSpace, output_path: Path) -> None:
    palette = ("#d9e9f4", "#91bad2", "#70a878", "#c77c75")
    fig, ax = plt.subplots(figsize=(5.2, 4.5), dpi=240)
    for niche_index, (name, color) in enumerate(zip(NICHE_NAMES, palette)):
        subset = simulation.meta[simulation.meta["niche"].astype(int) == niche_index]
        ax.scatter(
            subset["col"],
            -subset["row"],
            s=5.2,
            color=color,
            edgecolors="none",
            alpha=0.9,
            rasterized=True,
            label=name,
        )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.15),
        frameon=False,
        markerscale=2.2,
        fontsize=11,
        ncols=2,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _target_crop(target_path: Path) -> np.ndarray:
    image = np.asarray(Image.open(target_path).convert("RGB"))
    height, width = image.shape[:2]
    # Fractions reproduce the plot-only crop of the supplied 496 x 428 image
    # while remaining robust to a proportionally resized screenshot.
    y0, y1 = int(0.22 * height), height
    x0, x1 = int(0.25 * width), width
    crop = image[y0:y1, x0:x1].copy()
    # The legend overlaps the left edge of the point-cloud bounding box in the
    # supplied screenshot. Remove only dark grayscale glyphs, retaining the
    # colored (including pale-blue) point marks.
    chroma = np.ptp(crop.astype(int), axis=2)
    brightness = crop.mean(axis=2)
    crop[(chroma < 7) & (brightness < 248)] = 255
    return crop


def plot_comparison(
    simulation: ss.SimSpace,
    output_path: Path,
    target_path: Path,
    diagnostics: dict[str, float],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), dpi=240)
    axes[0].imshow(_target_crop(target_path))
    axes[0].set_title("Visual target (scHolography Fig. 3A)", fontsize=10)
    axes[0].axis("off")
    _scatter_simulation(axes[1], simulation, point_size=4.8)
    axes[1].set_title(
        f"Reference-free SimSpace (n={diagnostics['n_cells']:,})",
        fontsize=10,
    )
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=5,
        frameon=False,
        fontsize=8,
        markerscale=2.0,
    )
    fig.subplots_adjust(bottom=0.18, wspace=0.04)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def candidate_configs(base: CortexConfig, round_name: str) -> list[tuple[str, CortexConfig]]:
    if round_name == "geometry":
        return [
            ("G1 compact", replace(base, center_x=1.10, outer_radius_slope=0.15)),
            ("G2 selected", base),
            ("G3 broad", replace(base, outer_radius_at_top=1.17, outer_radius_slope=0.20)),
            ("G4 steep", replace(base, angle_min_deg=-166.0, center_y=1.22)),
            ("G5 tapered", replace(base, outer_radius_at_top=1.11, outer_radius_slope=0.23)),
            ("G6 long arc", replace(base, angle_min_deg=-172.0, center_x=1.18)),
        ]
    if round_name == "mrf":
        return [
            ("M1 phi 2.2", replace(base, phi=2.20, num_iterations=1, simulation_seed=11)),
            ("M2 phi 2.6", base),
            ("M3 phi 3.1", replace(base, phi=3.10, num_iterations=1, simulation_seed=41)),
            ("M4 phi 3.6", replace(base, phi=3.60, num_iterations=1, simulation_seed=53)),
            (
                "M5 local 2x",
                replace(
                    base,
                    neighbor_distance=1,
                    phi=2.70,
                    num_iterations=2,
                    simulation_seed=67,
                ),
            ),
            (
                "M6 broad",
                replace(
                    base,
                    neighbor_distance=3,
                    phi=2.70,
                    num_iterations=1,
                    simulation_seed=79,
                ),
            ),
        ]
    if round_name == "seed":
        seeds = (7, 17, 23, 41, 73, 101)
        return [
            (
                f"S{seed}",
                replace(base, simulation_seed=seed, initial_seed=1000 + seed),
            )
            for seed in seeds
        ]
    raise ValueError(f"Unknown sweep round: {round_name}")


def run_sweep(
    base: CortexConfig,
    round_name: str,
    output_dir: Path,
) -> pd.DataFrame:
    candidates = candidate_configs(base, round_name)
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.0), dpi=190)
    rows: list[dict[str, float | str]] = []
    for ax, (name, config) in zip(axes.ravel(), candidates):
        simulation, _, _ = simulate_cortex(config)
        _, diagnostics = summarize_simulation(simulation)
        _scatter_simulation(ax, simulation, point_size=2.8)
        ax.set_title(
            f"{name}\nscore={diagnostics['tuning_score']:.3f}, "
            f"n={diagnostics['n_cells']:,}",
            fontsize=9,
        )
        rows.append({"candidate": name, **asdict(config), **diagnostics})
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=5,
        frameon=False,
        fontsize=8,
        markerscale=2.0,
    )
    fig.subplots_adjust(bottom=0.12, hspace=0.18, wspace=0.05)
    fig.savefig(
        output_dir / f"candidate_sweep_{round_name}.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)
    results = pd.DataFrame(rows).sort_values("tuning_score")
    results.to_csv(output_dir / f"candidate_sweep_{round_name}.tsv", sep="\t", index=False)
    return results


def write_parameters(
    config: CortexConfig,
    output_path: Path,
    niche_mixture: np.ndarray = DEFAULT_NICHE_MIXTURE,
) -> None:
    theta = build_theta_matrices(
        niche_mixture,
        bias_strength=config.theta_bias_strength,
        depth_penalty=config.theta_depth_penalty,
        rare_state_cutoff=config.rare_state_cutoff,
        rare_state_penalty=config.rare_state_penalty,
    )
    payload = {
        "mode": "reference-free",
        "spatial_reference_used_for_generation": False,
        "n_group": len(NICHE_NAMES),
        "n_state": len(CELL_TYPES),
        "config": asdict(config),
        "cell_types": list(CELL_TYPES),
        "cell_colors": list(CELL_COLORS),
        "niche_names": list(NICHE_NAMES),
        "niche_cell_type_probabilities": niche_mixture.tolist(),
        "theta_matrices": [matrix.tolist() for matrix in theta],
        "visual_tuning_guides": {
            "approximate_composition": VISUAL_COMPOSITION_GUIDE.tolist(),
            "expected_median_cortical_depth": EXPECTED_DEPTH.tolist(),
            "formal_validation_metric": False,
        },
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
    )
    parser.add_argument(
        "--target-image",
        type=Path,
        help="Optional screenshot used only for the side-by-side visual panel.",
    )
    parser.add_argument(
        "--sweep",
        choices=("geometry", "mrf", "seed"),
        action="append",
        default=[],
        help="Generate one or more six-candidate tuning panels.",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> None:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config = CortexConfig()

    for round_name in args.sweep:
        results = run_sweep(config, round_name, args.output_dir)
        print(f"\n{round_name} sweep")
        print(results[["candidate", "tuning_score", "n_cells"]].to_string(index=False))

    simulation, _, _ = simulate_cortex(config)
    summary, diagnostics = summarize_simulation(simulation)
    simulation.meta.to_csv(args.output_dir / "simulated_cortex_cells.csv", index=False)
    summary.to_csv(args.output_dir / "cell_type_summary.tsv", sep="\t", index=False)
    (args.output_dir / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )
    write_parameters(config, args.output_dir / "cortex_parameters.json")
    plot_single(simulation, args.output_dir / "simulated_cortex.png")
    plot_niches(simulation, args.output_dir / "simulated_cortex_niches.png")
    if args.target_image:
        plot_comparison(
            simulation,
            args.output_dir / "target_vs_simulation.png",
            args.target_image,
            diagnostics,
        )
    print("\nfinal configuration")
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
