"""Supplementary Figure 5: Gibbs-sweep sensitivity analysis for R1-8.

This analysis varies only the phenotype-level ``num_iterations`` parameter in
the reference-free configuration used for Supplementary Figure 4. It
characterizes the sensitivity of declared spatial summaries; it is not a
formal diagnostic or proof of full-field equilibrium.
"""

from __future__ import annotations

import json
import math
import os
import platform
import subprocess
import sys
from copy import deepcopy
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/simspace-r1-8-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/simspace-r1-8-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

import simspace as ss


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "Panel_A_B_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig5"

SWEEP_COUNTS = (1, 2, 3, 4, 6, 8, 10)
DEFAULT_SWEEPS = 4
SEEDS = tuple(range(8))
SHAPE = (100, 100)
N_NICHES = 2
N_STATES = 9
NICHE_SWEEPS = 6
NEIGHBORHOOD_DISTANCE = 3
PERTURBATION_SD = 0.2
MORAN_K = 5
GEARY_K = 20
INTERACTION_K = 50
SIMSPACE_EXPECTED_VERSION = "0.3.4"
SIMSPACE_SOURCE_COMMIT = "de0a4c002e4ae733e354e3e180ab69b381ad994a"

METRIC_LABELS = {
    "moran_i": "Moran's I",
    "geary_c": "Geary's C",
    "interaction_distance": "Interaction distance",
}

METRIC_TITLES = {
    "moran_i": "Moran's I Distribution",
    "geary_c": "Geary's C Distribution",
    "interaction_distance": "Cell Type Interaction Distribution",
}


def package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "not installed as a distribution"


def adjacent_simspace_commit() -> str:
    source_checkout = REPO_ROOT.parent / "SimSpace"
    if not (source_checkout / ".git").exists():
        return "adjacent source checkout not found"
    result = subprocess.run(
        ["git", "-C", str(source_checkout), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def knn_edges(coordinates: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    effective_k = min(k, len(coordinates) - 1)
    neighbors = cKDTree(coordinates).query(
        coordinates, k=effective_k + 1
    )[1]
    if neighbors.ndim == 1:
        neighbors = neighbors[:, None]
    neighbors = neighbors[:, 1:]
    source = np.repeat(np.arange(len(coordinates)), neighbors.shape[1])
    target = neighbors.reshape(-1)
    return source, target


def moran_i(values: np.ndarray, edges: tuple[np.ndarray, np.ndarray]) -> float:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    denominator = float(centered @ centered)
    if denominator == 0:
        return np.nan
    source, target = edges
    return float(
        len(values)
        / len(source)
        * np.sum(centered[source] * centered[target])
        / denominator
    )


def geary_c(values: np.ndarray, edges: tuple[np.ndarray, np.ndarray]) -> float:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    denominator = float(centered @ centered)
    if denominator == 0:
        return np.nan
    source, target = edges
    return float(
        (len(values) - 1)
        / (2 * len(source))
        * np.sum((values[source] - values[target]) ** 2)
        / denominator
    )


def parameter_matrices(
    parameters: dict[str, object],
) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    niche_theta = np.zeros((N_NICHES, N_NICHES))
    niche_theta[np.triu_indices(N_NICHES, 1)] = parameters["niche_theta"]
    niche_theta = niche_theta + niche_theta.T
    np.fill_diagonal(niche_theta, 1)

    theta_list: list[np.ndarray] = []
    for niche in range(N_NICHES):
        theta = np.zeros((N_STATES, N_STATES))
        theta[np.triu_indices(N_STATES, 1)] = parameters["theta_list"][niche]
        theta = theta + theta.T
        np.fill_diagonal(theta, 1)
        theta_list.append(theta)

    density = np.asarray(parameters["density_replicates"], dtype=float).copy()
    density[density < 0] = 0
    return niche_theta, theta_list, density


def simulate_checkpoints(seed: int) -> dict[int, object]:
    """Run one Gibbs trajectory and materialize the declared sweep checkpoints."""
    parameters = ss.util.generate_random_parameters(
        n_group=N_NICHES,
        n_state=N_STATES,
        seed=seed,
    )
    niche_theta, theta_list, density = parameter_matrices(parameters)
    simulation = ss.SimSpace(
        shape=SHAPE,
        num_states=N_STATES,
        num_iterations=max(SWEEP_COUNTS),
        theta=theta_list,
        phi=float(parameters["phi_replicates"]),
        neighborhood=ss.spatial.generate_offsets(
            NEIGHBORHOOD_DISTANCE, "manhattan"
        ),
        random_seed=seed,
    )
    simulation.initialize()
    simulation.create_niche(
        num_niches=N_NICHES,
        n_iter=NICHE_SWEEPS,
        theta_niche=niche_theta,
    )

    checkpoints: dict[int, object] = {}
    np.random.seed(seed)
    for sweep in range(1, max(SWEEP_COUNTS) + 1):
        for row in np.random.permutation(SHAPE[0])[
            : math.ceil(SHAPE[0] * simulation.rho)
        ]:
            for col in np.random.permutation(SHAPE[1])[
                : math.ceil(SHAPE[1] * simulation.rho)
            ]:
                if simulation.grid[row, col] < 0:
                    continue
                neighbors = simulation.get_custom_neighbors(
                    row, col, neighborhood=simulation.neighborhood
                )
                if not neighbors:
                    continue
                niche_class = simulation.niche[row, col]
                probabilities = simulation._conditional_probability(
                    simulation.grid,
                    simulation.theta[niche_class],
                    neighbors,
                )
                simulation.grid[row, col] = np.random.choice(
                    range(N_STATES), p=probabilities
                )

        if sweep in SWEEP_COUNTS:
            rng_state = np.random.get_state()
            checkpoint = deepcopy(simulation)
            checkpoint.density_sampler(density)
            checkpoint.perturbation(step=PERTURBATION_SD)
            checkpoints[sweep] = checkpoint
            np.random.set_state(rng_state)

    return checkpoints


def measure_simulation(simulation, seed: int, sweeps: int) -> pd.DataFrame:
    coordinates = simulation.meta[["row", "col"]].to_numpy(dtype=float)
    labels = simulation.meta["state"].astype(int).to_numpy()
    moran_edges = knn_edges(coordinates, MORAN_K)
    geary_edges = knn_edges(coordinates, GEARY_K)
    rows: list[dict[str, object]] = []

    for state in range(N_STATES):
        indicator = (labels == state).astype(float)
        rows.extend(
            [
                {
                    "seed": seed,
                    "num_iterations": sweeps,
                    "metric": "moran_i",
                    "component": f"state_{state}",
                    "value": moran_i(indicator, moran_edges),
                },
                {
                    "seed": seed,
                    "num_iterations": sweeps,
                    "metric": "geary_c",
                    "component": f"state_{state}",
                    "value": geary_c(indicator, geary_edges),
                },
            ]
        )

    states = list(range(N_STATES))
    interaction = ss.spatial.calculate_interaction_score(
        pd.Series(labels),
        coordinates,
        states,
        k=INTERACTION_K,
    )
    for state_a in states:
        for state_b in states:
            rows.append(
                {
                    "seed": seed,
                    "num_iterations": sweeps,
                    "metric": "interaction_distance",
                    "component": f"state_{state_a}:state_{state_b}",
                    "value": float(interaction.loc[state_a, state_b]),
                }
            )

    return pd.DataFrame(rows)


def agreement_with_default(raw_metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for seed in SEEDS:
        for metric in METRIC_LABELS:
            subset = raw_metrics[
                (raw_metrics["seed"] == seed)
                & (raw_metrics["metric"] == metric)
            ]
            pivot = subset.pivot(
                index="component", columns="num_iterations", values="value"
            )
            reference = pivot[DEFAULT_SWEEPS]
            for sweeps in SWEEP_COUNTS:
                pair = pd.concat([reference, pivot[sweeps]], axis=1).dropna()
                pair.columns = ["default", "comparison"]
                correlation = (
                    1.0
                    if sweeps == DEFAULT_SWEEPS
                    else pearsonr(pair["default"], pair["comparison"]).statistic
                )
                rows.append(
                    {
                        "seed": seed,
                        "metric": metric,
                        "num_iterations": sweeps,
                        "pearson_r_vs_4": float(correlation),
                        "rmse_vs_4": float(
                            np.sqrt(
                                np.mean(
                                    (
                                        pair["comparison"].to_numpy()
                                        - pair["default"].to_numpy()
                                    )
                                    ** 2
                                )
                            )
                        ),
                        "n_components": len(pair),
                    }
                )
    return pd.DataFrame(rows)


def summarize_distributions(raw_metrics: pd.DataFrame) -> pd.DataFrame:
    return (
        raw_metrics.groupby(["metric", "num_iterations"])["value"]
        .agg(
            n="count",
            mean="mean",
            sd="std",
            median="median",
            q25=lambda values: values.quantile(0.25),
            q75=lambda values: values.quantile(0.75),
        )
        .reset_index()
    )


def summarize_agreement(agreement: pd.DataFrame) -> pd.DataFrame:
    return (
        agreement.groupby(["metric", "num_iterations"])
        .agg(
            median_pearson_r_vs_4=("pearson_r_vs_4", "median"),
            min_pearson_r_vs_4=("pearson_r_vs_4", "min"),
            max_pearson_r_vs_4=("pearson_r_vs_4", "max"),
            median_rmse_vs_4=("rmse_vs_4", "median"),
        )
        .reset_index()
    )


def plot_figure(
    raw_metrics: pd.DataFrame,
    representative: dict[int, pd.DataFrame],
) -> plt.Figure:
    sns.set_theme(style="white", context="paper")
    figure = plt.figure(figsize=(12.0, 6.7), dpi=300)
    grid = figure.add_gridspec(2, 1, height_ratios=[1.0, 1.0])
    summary_grid = grid[0, 0].subgridspec(1, 3, wspace=0.30)
    layout_grid = grid[1, 0].subgridspec(1, 4, wspace=0.25)
    palette = sns.color_palette("viridis", len(SWEEP_COUNTS))

    for index, metric in enumerate(METRIC_LABELS):
        axis = figure.add_subplot(summary_grid[0, index])
        subset = raw_metrics[raw_metrics["metric"] == metric]
        sns.violinplot(
            data=subset,
            x="num_iterations",
            y="value",
            order=SWEEP_COUNTS,
            hue="num_iterations",
            palette=palette,
            legend=False,
            linewidth=0.8,
            ax=axis,
        )
        axis.set_xlabel("Phenotype Gibbs sweeps")
        axis.set_ylabel(METRIC_LABELS[metric])
        axis.set_title(METRIC_TITLES[metric])
        if index == 0:
            axis.text(
                -0.12,
                1.08,
                "A",
                transform=axis.transAxes,
                fontsize=13,
                fontweight="bold",
                va="top",
            )

    state_palette = sns.color_palette("tab20", N_STATES)
    for index, sweeps in enumerate((1, 2, 4, 10)):
        axis = figure.add_subplot(layout_grid[0, index])
        meta = representative[sweeps]
        colors = [state_palette[int(state)] for state in meta["state"]]
        axis.scatter(
            meta["col"],
            meta["row"],
            c=colors,
            s=8,
            linewidths=0,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_title(
            f"{chr(ord('B') + index)}   {sweeps} sweep"
            f"{'s' if sweeps != 1 else ''}",
            loc="left",
            fontweight="bold",
        )
        axis.set_xlabel("X")
        axis.set_ylabel("Y")

    figure.tight_layout()
    return figure


def write_provenance() -> None:
    simspace_commit = adjacent_simspace_commit()
    if simspace_commit != SIMSPACE_SOURCE_COMMIT:
        raise RuntimeError(
            "Unexpected adjacent SimSpace source commit: "
            f"{simspace_commit}; expected {SIMSPACE_SOURCE_COMMIT}."
        )
    simspace_version = package_version("simspace")
    if simspace_version != SIMSPACE_EXPECTED_VERSION:
        raise RuntimeError(
            f"Unexpected SimSpace version: {simspace_version}; "
            f"expected {SIMSPACE_EXPECTED_VERSION}."
        )

    configuration = {
        "analysis": "R1-8 phenotype Gibbs-sweep sensitivity",
        "scope_note": (
            "Sensitivity of declared spatial summaries; not a formal "
            "diagnostic or proof of full-field equilibrium."
        ),
        "num_iterations": list(SWEEP_COUNTS),
        "default_num_iterations": DEFAULT_SWEEPS,
        "seeds": list(SEEDS),
        "shape": list(SHAPE),
        "num_niches": N_NICHES,
        "num_states": N_STATES,
        "fixed_niche_sweeps": NICHE_SWEEPS,
        "neighborhood_distance": NEIGHBORHOOD_DISTANCE,
        "neighborhood_method": "manhattan",
        "perturbation_sd": PERTURBATION_SD,
        "moran_knn": MORAN_K,
        "geary_knn": GEARY_K,
        "interaction_knn": INTERACTION_K,
        "simspace_version": simspace_version,
        "simspace_source_commit": simspace_commit,
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
            {"software": "matplotlib", "version": matplotlib.__version__},
            {"software": "seaborn", "version": sns.__version__},
            {"software": "platform", "version": platform.platform()},
            {"software": "python_executable", "version": sys.executable},
        ]
    )
    software.to_csv(DATA_DIR / "software_versions.tsv", sep="\t", index=False)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    write_provenance()

    metric_frames: list[pd.DataFrame] = []
    representative: dict[int, pd.DataFrame] = {}
    for seed in SEEDS:
        print(f"seed={seed}", flush=True)
        checkpoints = simulate_checkpoints(seed)
        for sweeps in SWEEP_COUNTS:
            simulation = checkpoints[sweeps]
            metric_frames.append(measure_simulation(simulation, seed, sweeps))
            if seed == 0 and sweeps in (1, 2, 4, 10):
                representative[sweeps] = simulation.meta.copy()

    raw_metrics = pd.concat(metric_frames, ignore_index=True)
    distribution_summary = summarize_distributions(raw_metrics)
    agreement = agreement_with_default(raw_metrics)
    agreement_summary = summarize_agreement(agreement)

    raw_metrics.to_csv(DATA_DIR / "spatial_summary_metrics.tsv", sep="\t", index=False)
    distribution_summary.to_csv(
        DATA_DIR / "distribution_summary.tsv", sep="\t", index=False
    )
    agreement.to_csv(
        DATA_DIR / "agreement_with_default.tsv", sep="\t", index=False
    )
    agreement_summary.to_csv(
        DATA_DIR / "agreement_summary.tsv", sep="\t", index=False
    )

    figure = plot_figure(raw_metrics, representative)
    output_path = SCRIPT_DIR / "SFig5.png"
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    figure.savefig(EXAMPLE_DIR / "SFig5.png", dpi=300, bbox_inches="tight")
    plt.close(figure)

    late = agreement_summary[
        agreement_summary["num_iterations"].isin((6, 8, 10))
    ]
    print("\nAgreement with the four-sweep default at later sweep counts:")
    print(late.to_string(index=False))


if __name__ == "__main__":
    main()
