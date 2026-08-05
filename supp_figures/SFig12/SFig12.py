"""Supplementary Figure 12: optional spatial-expression and observation layers.

This reproducible analysis addresses R2-M1 and R2-m5. It validates the
reference-free coordinate-conditioned log-mean option across predeclared
effect sizes and separately validates capture thinning, ambient background,
and mean-dependent dropout. The resulting conditional gene truth is also
structured for reuse in the related R2-M3 SVG benchmark.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/simspace-r2-expression-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/simspace-r2-expression-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

import simspace as ss


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "Panel_A_D_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig12"

SEEDS = tuple(range(20))
SHAPE = (50, 50)
N_STATES = 4
N_GENES = 500
MRF_SWEEPS = 3
MRF_PHI = 1.8
DENSITY_RETENTION = 0.65
COORDINATE_JITTER = 0.55
SPATIAL_STRENGTHS = (-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0)
NONZERO_SPATIAL_STRENGTHS = tuple(
    strength for strength in SPATIAL_STRENGTHS if strength != 0.0
)
N_EFFECT_GENES = 96
N_NULL_VALIDATION_GENES = 24
SPATIAL_DIRECTION = np.asarray((1.0, 1.0), dtype=float) / np.sqrt(2.0)
LINEAR_DIRECTIONS = {
    "Horizontal": np.asarray((1.0, 0.0), dtype=float),
    "Vertical": np.asarray((0.0, 1.0), dtype=float),
    "Diagonal": SPATIAL_DIRECTION,
    "Anti-diagonal": np.asarray((1.0, -1.0), dtype=float) / np.sqrt(2.0),
}
SPATIAL_BASIS_ORDER = ("linear", "radial", "hotspot", "structure_distance")
SPATIAL_BASIS_LABELS = {
    "linear": "Linear",
    "radial": "Radial",
    "hotspot": "Hotspot",
    "structure_distance": "Structure distance",
}
SPATIAL_BASIS_COLORS = {
    "linear": "#4C78A8",
    "radial": "#F58518",
    "hotspot": "#54A24B",
    "structure_distance": "#B279A2",
}
SPATIAL_BASIS_MARKERS = {
    "linear": "o",
    "radial": "s",
    "hotspot": "^",
    "structure_distance": "D",
}
CAPTURE_LEVELS = (1.0, 0.85, 0.7, 0.55, 0.4)
AMBIENT_LEVELS = (0.0, 5.0, 20.0, 50.0, 100.0)
DROPOUT_LEVELS = (
    ("Off", None),
    ("Very low", {"mode": "mean_dependent", "intercept": -0.5, "slope": -1.0}),
    ("Low", {"mode": "mean_dependent", "intercept": 0.5, "slope": -1.0}),
    ("Moderate", {"mode": "mean_dependent", "intercept": 1.5, "slope": -1.0}),
    ("High", {"mode": "mean_dependent", "intercept": 2.5, "slope": -1.0}),
)
COMBINED_TECHNICAL_NOISE = {
    "capture_efficiency": 0.7,
    "ambient_rate": 50.0,
}
COMBINED_DROPOUT = {
    "mode": "mean_dependent",
    "intercept": 1.5,
    "slope": -1.0,
}
SIMSPACE_EXPECTED_VERSION = "0.3.4"
SIMSPACE_SOURCE_COMMIT = "de0a4c002e4ae733e354e3e180ab69b381ad994a"
REPRESENTATIVE_SEED = 0
REPRESENTATIVE_GENE_ID = 95
REPRESENTATIVE_GENE = f"Gene_{REPRESENTATIVE_GENE_ID}"
REPRESENTATIVE_GENE_SELECTION = (
    "Among seed-0 linear genes with configured strength +1, select the gene "
    "with the largest minimum phenotype-specific baseline mean; break ties by "
    "the median phenotype-specific baseline mean and then GeneID."
)

CONDITION_ORDER = ("Clean", "Capture 0.4", "Ambient 100", "Dropout high", "Combined")
CONDITION_COLORS = {
    "Clean": "#4C78A8",
    "Capture 0.4": "#59A14F",
    "Ambient 100": "#F28E2B",
    "Dropout high": "#B279A2",
    "Combined": "#E15759",
}
STATE_COLORS = ("#4C78A8", "#F28E2B", "#59A14F", "#B279A2")


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


def adjacent_simspace_source_dirty() -> bool | str:
    """Report uncommitted package-source changes, excluding unrelated files."""
    source_checkout = REPO_ROOT.parent / "SimSpace"
    if not (source_checkout / ".git").exists():
        return "adjacent source checkout not found"
    result = subprocess.run(
        [
            "git",
            "-C",
            str(source_checkout),
            "status",
            "--short",
            "--",
            "simspace",
            "pyproject.toml",
            "setup.py",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return bool(result.stdout.strip())


def fixed_theta() -> list[np.ndarray]:
    """Return one fixed phenotype-affinity matrix for the single niche."""
    theta = np.full((N_STATES, N_STATES), 0.08, dtype=float)
    np.fill_diagonal(theta, 1.0)
    return [theta]


def build_layout(seed: int) -> pd.DataFrame:
    """Generate one reference-free spatial label layout."""
    simulation = ss.SimSpace(
        shape=SHAPE,
        num_states=N_STATES,
        num_iterations=MRF_SWEEPS,
        theta=fixed_theta(),
        phi=MRF_PHI,
        neighborhood=ss.spatial.generate_offsets(1, "manhattan"),
        random_seed=1000 + seed,
    )
    simulation.initialize()
    simulation.gibbs_sampler()
    simulation.density_sampler(float(DENSITY_RETENTION))
    simulation.perturbation(step=COORDINATE_JITTER)
    counts = simulation.meta["state"].astype(int).value_counts()
    if len(counts) != N_STATES or counts.min() < 50:
        raise RuntimeError(f"Seed {seed} did not retain all four cell types adequately.")
    return simulation.meta[["state", "row", "col", "niche"]].copy()


def expression_simulation(meta: pd.DataFrame, seed: int) -> ss.SimSpace:
    simulation = ss.SimSpace(
        shape=SHAPE,
        num_states=N_STATES,
        random_seed=10000 + seed,
    )
    simulation.meta = meta.copy()
    return simulation


def validation_effect_design(basis_name: str) -> pd.DataFrame:
    """Declare spatial-effect and null-control genes for one basis family."""
    if basis_name not in SPATIAL_BASIS_ORDER:
        raise ValueError(f"Unknown spatial basis: {basis_name}")

    rows: list[dict[str, object]] = []
    if basis_name == "linear":
        combinations = [
            (strength, direction_name, replicate)
            for strength in NONZERO_SPATIAL_STRENGTHS
            for direction_name in LINEAR_DIRECTIONS
            for replicate in range(3)
        ]
        representative = (1.0, "Diagonal", 0)
        combinations.remove(representative)
        combinations.insert(0, representative)
        for gene_id, (strength, direction_name, replicate) in enumerate(combinations):
            direction = LINEAR_DIRECTIONS[direction_name]
            rows.append(
                {
                    "GeneID": gene_id,
                    "Gene": f"Gene_{gene_id}",
                    "configured_strength": float(strength),
                    "direction": direction_name,
                    "direction_row": float(direction[0]),
                    "direction_col": float(direction[1]),
                    "replicate": replicate,
                    "is_direct": True,
                }
            )
    else:
        combinations = [
            (strength, replicate)
            for strength in NONZERO_SPATIAL_STRENGTHS
            for replicate in range(12)
        ]
        for gene_id, (strength, replicate) in enumerate(combinations):
            rows.append(
                {
                    "GeneID": gene_id,
                    "Gene": f"Gene_{gene_id}",
                    "configured_strength": float(strength),
                    "direction": "Not applicable",
                    "direction_row": np.nan,
                    "direction_col": np.nan,
                    "replicate": replicate,
                    "is_direct": True,
                }
            )

    if len(rows) != N_EFFECT_GENES:
        raise RuntimeError("The declared spatial-effect design has the wrong size.")

    direction_names = tuple(LINEAR_DIRECTIONS)
    for null_index in range(N_NULL_VALIDATION_GENES):
        direction_name = (
            direction_names[null_index % len(direction_names)]
            if basis_name == "linear"
            else "Not applicable"
        )
        direction = (
            LINEAR_DIRECTIONS[direction_name]
            if basis_name == "linear"
            else np.asarray((np.nan, np.nan))
        )
        gene_id = N_EFFECT_GENES + null_index
        rows.append(
            {
                "GeneID": gene_id,
                "Gene": f"Gene_{gene_id}",
                "configured_strength": 0.0,
                "direction": direction_name,
                "direction_row": float(direction[0]),
                "direction_col": float(direction[1]),
                "replicate": null_index,
                "is_direct": False,
            }
        )
    return pd.DataFrame(rows)


def direct_spatial_config(
    basis_name: str,
    *,
    zero_coefficients: bool = False,
) -> dict[str, object]:
    """Build one supported basis configuration with heterogeneous gene effects."""
    design = validation_effect_design(basis_name)
    direct = design[design["is_direct"]].copy()
    coefficients: dict[str, list[float]] = {}
    for row in direct.itertuples(index=False):
        strength = 0.0 if zero_coefficients else float(row.configured_strength)
        if basis_name == "linear":
            coefficients[row.Gene] = [
                strength * float(row.direction_row),
                strength * float(row.direction_col),
            ]
        else:
            coefficients[row.Gene] = [strength]

    config: dict[str, object] = {
        "genes": direct["Gene"].tolist(),
        "basis": basis_name,
        "coefficients": coefficients,
        "coordinate_scaling": "unit_range",
        "center_basis": True,
    }
    if basis_name == "radial":
        config["center"] = [0.5, 0.5]
    elif basis_name == "hotspot":
        config["center"] = [0.30, 0.70]
        config["bandwidth"] = 0.18
    elif basis_name == "structure_distance":
        config["structure_points"] = [
            [0.15, 0.20],
            [0.325, 0.35],
            [0.50, 0.50],
            [0.675, 0.65],
            [0.85, 0.80],
        ]
    return config


def simulate_expression(
    meta: pd.DataFrame,
    seed: int,
    *,
    spatial_basis: str | None,
    zero_spatial: bool = False,
    technical_noise: dict[str, object] | None = None,
    dropout: dict[str, object] | None = None,
) -> ss.SimSpace:
    simulation = expression_simulation(meta, seed)
    simulation.create_omics(
        n_genes=N_GENES,
        bg_ratio=1.0,
        bg_param=(6.0, 2.0),
        marker_param=(6.0, 2.0),
        lr_ratio=0.0,
        spatial=False,
        direct_spatial_effects=(
            None
            if spatial_basis is None
            else direct_spatial_config(spatial_basis, zero_coefficients=zero_spatial)
        ),
        technical_noise=technical_noise,
        dropout=dropout,
        store_intermediate=True,
    )
    return simulation


def spatial_score(simulation: ss.SimSpace) -> np.ndarray:
    """Return the representative diagonal linear score used for Gene 0."""
    design = simulation.omics_spatial_design
    basis = design[["row_linear", "col_linear"]].to_numpy(dtype=float)
    return basis @ SPATIAL_DIRECTION


def validation_score(simulation: ss.SimSpace, effect_row: object, basis_name: str) -> np.ndarray:
    """Return the centered scalar score associated with one declared coefficient."""
    design = simulation.omics_spatial_design
    if basis_name == "linear":
        basis = design[["row_linear", "col_linear"]].to_numpy(dtype=float)
        direction = np.asarray(
            [float(effect_row.direction_row), float(effect_row.direction_col)],
            dtype=float,
        )
        return basis @ direction
    basis_column = {
        "radial": "radial_distance",
        "hotspot": "hotspot_basis",
        "structure_distance": "structure_distance",
    }[basis_name]
    return design[basis_column].to_numpy(dtype=float)


def poisson_offset_slope(
    counts: np.ndarray,
    baseline_mean: np.ndarray,
    score: np.ndarray,
) -> float:
    """Estimate the direct log-mean coefficient with a known baseline offset."""
    counts = np.asarray(counts, dtype=float)
    baseline_mean = np.asarray(baseline_mean, dtype=float)
    score = np.asarray(score, dtype=float)
    if np.any(baseline_mean <= 0):
        raise ValueError("The Poisson offset requires positive baseline means.")
    design = np.column_stack((np.ones(len(score)), score))
    offset = np.log(baseline_mean)
    coefficient = np.zeros(2, dtype=float)
    for _ in range(100):
        linear_predictor = np.clip(offset + design @ coefficient, -30.0, 30.0)
        fitted_mean = np.exp(linear_predictor)
        score_vector = design.T @ (counts - fitted_mean)
        information = design.T @ (fitted_mean[:, None] * design)
        step = np.linalg.solve(information + np.eye(2) * 1e-10, score_vector)
        coefficient += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return float(coefficient[1])


def gene_latent_and_observed_slopes(
    simulation: ss.SimSpace,
    gene_name: str,
    score: np.ndarray,
) -> tuple[float, float, pd.DataFrame, np.ndarray]:
    baseline = ss.omics.buildOmicsMean(simulation.gene_meta, simulation.meta)
    latent_log_ratio = np.log(
        simulation.omics_latent_mean[gene_name].to_numpy(dtype=float)
        / baseline[gene_name].to_numpy(dtype=float)
    )
    latent_design = np.column_stack((np.ones(len(score)), score))
    latent_slope = float(np.linalg.lstsq(latent_design, latent_log_ratio, rcond=None)[0][1])
    observed_slope = poisson_offset_slope(
        simulation.omics[gene_name].to_numpy(dtype=float),
        baseline[gene_name].to_numpy(dtype=float),
        score,
    )
    return latent_slope, observed_slope, baseline, score


def latent_and_observed_slopes(
    simulation: ss.SimSpace,
) -> tuple[float, float, pd.DataFrame, np.ndarray]:
    """Recover the representative Gene 0 diagonal effect."""
    score = spatial_score(simulation)
    return gene_latent_and_observed_slopes(simulation, "Gene_0", score)


def expression_metrics(
    simulation: ss.SimSpace,
    condition: str,
    seed: int,
) -> dict[str, object]:
    latent_slope, observed_slope, _, _ = latent_and_observed_slopes(simulation)
    counts = simulation.omics.to_numpy(dtype=float)
    latent = simulation.omics_latent_mean.to_numpy(dtype=float)
    gene_means = counts.mean(axis=0)
    gene_variances = counts.var(axis=0, ddof=1)
    mean_variance_rho = float(spearmanr(gene_means, gene_variances).statistic)
    latent_recovery = float(
        pearsonr(np.log1p(counts).ravel(), np.log1p(latent).ravel()).statistic
    )
    return {
        "seed": seed,
        "condition": condition,
        "median_library_size": float(np.median(counts.sum(axis=1))),
        "zero_fraction": float(np.mean(counts == 0)),
        "mean_variance_spearman": mean_variance_rho,
        "latent_recovery_pearson": latent_recovery,
        "latent_spatial_slope": latent_slope,
        "observed_spatial_slope": observed_slope,
    }


def gene_truth_rows(simulation: ss.SimSpace, seed: int) -> list[dict[str, object]]:
    effect_design = validation_effect_design("linear").set_index("GeneID")
    baseline = ss.omics.buildOmicsMean(simulation.gene_meta, simulation.meta)
    rows: list[dict[str, object]] = []
    type_columns = [f"Type_{state}" for state in range(N_STATES)]
    for _, gene in simulation.gene_meta.iterrows():
        gene_id = int(gene["GeneID"])
        is_direct = gene_id in effect_design.index and bool(effect_design.loc[gene_id, "is_direct"])
        if is_direct:
            effect = effect_design.loc[gene_id]
            coefficient_row = float(effect["configured_strength"] * effect["direction_row"])
            coefficient_col = float(effect["configured_strength"] * effect["direction_col"])
            log_effect = np.log(
                simulation.omics_latent_mean[f"Gene_{gene_id}"].to_numpy(dtype=float)
                / baseline[f"Gene_{gene_id}"].to_numpy(dtype=float)
            )
        else:
            effect = None
            coefficient_row = 0.0
            coefficient_col = 0.0
            log_effect = np.zeros(len(simulation.meta), dtype=float)
        row: dict[str, object] = {
            "seed": seed,
            "GeneID": gene_id,
            "Gene": f"Gene_{gene_id}",
            "Marker": int(gene["Marker"]),
            "truth_scope": (
                "coordinate_conditioned_within_cell_type" if is_direct else "cell_type_only"
            ),
            "conditional_spatial_truth": is_direct,
            "spatial_basis": "linear" if is_direct else "none",
            "configured_strength": float(effect["configured_strength"]) if is_direct else 0.0,
            "direction": str(effect["direction"]) if is_direct else "none",
            "coefficient_row": coefficient_row,
            "coefficient_col": coefficient_col,
            "min_log_effect": float(log_effect.min()),
            "max_log_effect": float(log_effect.max()),
        }
        for column in type_columns:
            row[column] = float(gene[column])
        rows.append(row)
    return rows


def select_representative_gene(simulation: ss.SimSpace) -> int:
    """Select a visible example without using its realized count pattern."""
    design = validation_effect_design("linear")
    candidates = design[(design["is_direct"]) & (design["configured_strength"] == 1.0)]
    type_columns = [f"Type_{state}" for state in range(N_STATES)]
    gene_meta = simulation.gene_meta.set_index("GeneID")
    scores: list[dict[str, float | int]] = []
    for effect in candidates.itertuples(index=False):
        baseline_means = gene_meta.loc[int(effect.GeneID), type_columns].to_numpy(dtype=float)
        scores.append(
            {
                "GeneID": int(effect.GeneID),
                "minimum_baseline_mean": float(baseline_means.min()),
                "median_baseline_mean": float(np.median(baseline_means)),
            }
        )
    ranking = pd.DataFrame(scores).sort_values(
        ["minimum_baseline_mean", "median_baseline_mean", "GeneID"],
        ascending=[False, False, True],
    )
    return int(ranking.iloc[0]["GeneID"])


def representative_map(simulation: ss.SimSpace) -> pd.DataFrame:
    selected_gene_id = select_representative_gene(simulation)
    if selected_gene_id != REPRESENTATIVE_GENE_ID:
        raise RuntimeError(
            "The deterministic representative-gene rule no longer selects "
            f"{REPRESENTATIVE_GENE}; selected Gene_{selected_gene_id} instead."
        )
    baseline = ss.omics.buildOmicsMean(simulation.gene_meta, simulation.meta)
    latent = simulation.omics_latent_mean[REPRESENTATIVE_GENE].to_numpy(dtype=float)
    base = baseline[REPRESENTATIVE_GENE].to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "state": simulation.meta["state"].astype(int).to_numpy(),
            "row": simulation.meta["row"].to_numpy(dtype=float),
            "col": simulation.meta["col"].to_numpy(dtype=float),
            "representative_gene": REPRESENTATIVE_GENE,
            "direct_log_fold_truth": np.log(latent / base),
            "latent_mean": latent,
            "observed_count": simulation.omics[REPRESENTATIVE_GENE].to_numpy(dtype=float),
        }
    )


def summarize_values(
    frame: pd.DataFrame,
    group_columns: list[str],
    value_columns: list[str],
    section: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    grouped = frame.groupby(group_columns, dropna=False, sort=False)
    for group_key, group in grouped:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        group_label = "; ".join(
            f"{column}={value}" for column, value in zip(group_columns, group_key, strict=True)
        )
        for value_column in value_columns:
            values = group[value_column].to_numpy(dtype=float)
            rows.append(
                {
                    "section": section,
                    "group": group_label,
                    "metric": value_column,
                    "n_seeds": int(len(values)),
                    "median": float(np.median(values)),
                    "minimum": float(np.min(values)),
                    "maximum": float(np.max(values)),
                }
            )
    return rows


def acceptance_checks(
    effect: pd.DataFrame,
    calibration: pd.DataFrame,
    compatibility: pd.DataFrame,
) -> pd.DataFrame:
    seed_summary = (
        effect.groupby(["seed", "basis", "configured_strength"], as_index=False)[
            ["latent_strength", "observed_strength"]
        ]
        .median()
    )
    basis_medians = seed_summary.groupby(["basis", "configured_strength"])[
        "observed_strength"
    ].median()
    monotonic_differences = np.concatenate(
        [
            np.diff(
                basis_medians.xs(basis_name, level="basis")
                .sort_index()
                .to_numpy(dtype=float)
            )
            for basis_name in SPATIAL_BASIS_ORDER
        ]
    )
    configured_values = basis_medians.index.get_level_values("configured_strength").to_numpy(
        dtype=float
    )
    nonzero = seed_summary[seed_summary["configured_strength"] != 0].copy()
    sign_recovery = np.mean(
        np.sign(nonzero["configured_strength"].to_numpy())
        == np.sign(nonzero["observed_strength"].to_numpy())
    )
    null_medians = basis_medians.xs(0.0, level="configured_strength")
    capture = calibration[calibration["component"] == "capture"]
    ambient = calibration[calibration["component"] == "ambient"]
    dropout = calibration[calibration["component"] == "dropout"]

    checks = [
        {
            "check": "default_equals_disabled_components",
            "value": float(compatibility["default_equals_disabled"].mean()),
            "criterion": "equal to 1.0",
            "passed": bool(compatibility["default_equals_disabled"].all()),
        },
        {
            "check": "default_equals_zero_spatial_effect",
            "value": float(compatibility["default_equals_zero_spatial"].mean()),
            "criterion": "equal to 1.0",
            "passed": bool(compatibility["default_equals_zero_spatial"].all()),
        },
        {
            "check": "latent_spatial_coefficient_exact",
            "value": float(
                np.max(np.abs(effect["latent_strength"] - effect["configured_strength"]))
            ),
            "criterion": "maximum absolute error < 1e-10",
            "passed": bool(
                np.max(np.abs(effect["latent_strength"] - effect["configured_strength"]))
                < 1e-10
            ),
        },
        {
            "check": "observed_dose_response_monotonic",
            "value": float(np.min(monotonic_differences)),
            "criterion": "all within-basis adjacent median differences > 0",
            "passed": bool(np.all(monotonic_differences > 0)),
        },
        {
            "check": "observed_sign_recovery",
            "value": float(sign_recovery),
            "criterion": "at least 0.95 across seed-level gene summaries",
            "passed": bool(sign_recovery >= 0.95),
        },
        {
            "check": "observed_median_bias",
            "value": float(
                np.max(np.abs(configured_values - basis_medians.to_numpy(dtype=float)))
            ),
            "criterion": "maximum absolute median bias <= 0.15",
            "passed": bool(
                np.max(np.abs(configured_values - basis_medians.to_numpy(dtype=float)))
                <= 0.15
            ),
        },
        {
            "check": "observed_null_control_bias",
            "value": float(np.max(np.abs(null_medians.to_numpy(dtype=float)))),
            "criterion": "maximum absolute null-control median <= 0.10",
            "passed": bool(np.max(np.abs(null_medians.to_numpy(dtype=float))) <= 0.10),
        },
        {
            "check": "capture_calibration",
            "value": float(
                np.max(
                    np.abs(
                        capture.groupby("configured_value")["realized_value"].median()
                        - capture.groupby("configured_value")["configured_value"].median()
                    )
                )
            ),
            "criterion": "maximum median absolute error <= 0.03",
            "passed": bool(
                np.max(
                    np.abs(
                        capture.groupby("configured_value")["realized_value"].median()
                        - capture.groupby("configured_value")["configured_value"].median()
                    )
                )
                <= 0.03
            ),
        },
        {
            "check": "ambient_calibration",
            "value": float(
                np.max(
                    np.abs(
                        ambient.groupby("configured_value")["realized_value"].median()
                        - ambient.groupby("configured_value")["configured_value"].median()
                    )
                )
            ),
            "criterion": "maximum median absolute error <= 0.50 counts per cell",
            "passed": bool(
                np.max(
                    np.abs(
                        ambient.groupby("configured_value")["realized_value"].median()
                        - ambient.groupby("configured_value")["configured_value"].median()
                    )
                )
                <= 0.50
            ),
        },
        {
            "check": "dropout_calibration",
            "value": float(np.max(np.abs(dropout["realized_value"] - dropout["configured_value"]))),
            "criterion": "maximum seed-level absolute error <= 0.01",
            "passed": bool(
                np.max(np.abs(dropout["realized_value"] - dropout["configured_value"]))
                <= 0.01
            ),
        },
    ]
    return pd.DataFrame(checks)


def generate_data() -> dict[str, pd.DataFrame]:
    installed_version = package_version("simspace")
    if installed_version != SIMSPACE_EXPECTED_VERSION:
        raise RuntimeError(
            f"Supplementary Figure 12 requires SimSpace {SIMSPACE_EXPECTED_VERSION}; "
            f"found {installed_version}."
        )

    effect_rows: list[dict[str, object]] = []
    calibration_rows: list[dict[str, object]] = []
    observation_rows: list[dict[str, object]] = []
    compatibility_rows: list[dict[str, object]] = []
    truth_rows: list[dict[str, object]] = []
    representative: pd.DataFrame | None = None

    for seed in SEEDS:
        meta = build_layout(seed)
        baseline = simulate_expression(meta, seed, spatial_basis=None)
        disabled = simulate_expression(
            meta,
            seed,
            spatial_basis=None,
            technical_noise={},
            dropout={},
        )

        basis_simulations: dict[str, ss.SimSpace] = {}
        for basis_name in SPATIAL_BASIS_ORDER:
            simulation = simulate_expression(meta, seed, spatial_basis=basis_name)
            basis_simulations[basis_name] = simulation
            design = validation_effect_design(basis_name)
            baseline_mean = ss.omics.buildOmicsMean(simulation.gene_meta, simulation.meta)
            for effect in design.itertuples(index=False):
                score = validation_score(simulation, effect, basis_name)
                gene_name = str(effect.Gene)
                latent_log_ratio = np.log(
                    simulation.omics_latent_mean[gene_name].to_numpy(dtype=float)
                    / baseline_mean[gene_name].to_numpy(dtype=float)
                )
                latent_design = np.column_stack((np.ones(len(score)), score))
                latent_slope = float(
                    np.linalg.lstsq(latent_design, latent_log_ratio, rcond=None)[0][1]
                )
                observed_slope = poisson_offset_slope(
                    simulation.omics[gene_name].to_numpy(dtype=float),
                    baseline_mean[gene_name].to_numpy(dtype=float),
                    score,
                )
                effect_rows.append(
                    {
                        "seed": seed,
                        "basis": basis_name,
                        "basis_label": SPATIAL_BASIS_LABELS[basis_name],
                        "GeneID": int(effect.GeneID),
                        "Gene": gene_name,
                        "is_direct": bool(effect.is_direct),
                        "direction": str(effect.direction),
                        "configured_strength": float(effect.configured_strength),
                        "coefficient_row": (
                            float(effect.configured_strength * effect.direction_row)
                            if basis_name == "linear"
                            else np.nan
                        ),
                        "coefficient_col": (
                            float(effect.configured_strength * effect.direction_col)
                            if basis_name == "linear"
                            else np.nan
                        ),
                        "latent_strength": latent_slope,
                        "observed_strength": observed_slope,
                        "n_cells": len(meta),
                    }
                )

        zero_spatial = simulate_expression(
            meta,
            seed,
            spatial_basis="linear",
            zero_spatial=True,
        )
        compatibility_rows.append(
            {
                "seed": seed,
                "default_equals_disabled": bool(baseline.omics.equals(disabled.omics)),
                "default_equals_zero_spatial": bool(baseline.omics.equals(zero_spatial.omics)),
                "gene_metadata_equal_disabled": bool(
                    baseline.gene_meta.equals(disabled.gene_meta)
                ),
                "gene_metadata_equal_zero_spatial": bool(
                    baseline.gene_meta.equals(zero_spatial.gene_meta)
                ),
                "default_disabled_mismatched_entries": int(
                    np.count_nonzero(
                        baseline.omics.to_numpy() != disabled.omics.to_numpy()
                    )
                ),
                "default_zero_spatial_mismatched_entries": int(
                    np.count_nonzero(
                        baseline.omics.to_numpy() != zero_spatial.omics.to_numpy()
                    )
                ),
            }
        )

        clean = basis_simulations["linear"]
        truth_rows.extend(gene_truth_rows(clean, seed))
        condition_simulations: dict[str, ss.SimSpace] = {"Clean": clean}

        for efficiency in CAPTURE_LEVELS:
            if efficiency == 1.0:
                simulation = clean
            else:
                simulation = simulate_expression(
                    meta,
                    seed,
                    spatial_basis="linear",
                    technical_noise={
                        "capture_efficiency": efficiency,
                        "ambient_rate": 0.0,
                        "seed": 30000 + seed,
                    },
                )
            retained_fraction = float(
                simulation.omics_capture_counts.to_numpy(dtype=float).sum()
                / simulation.omics_baseline_counts.to_numpy(dtype=float).sum()
            )
            calibration_rows.append(
                {
                    "seed": seed,
                    "component": "capture",
                    "level": f"efficiency_{efficiency:g}",
                    "configured_value": efficiency,
                    "realized_value": retained_fraction,
                }
            )
            if efficiency == 0.4:
                condition_simulations["Capture 0.4"] = simulation

        for rate in AMBIENT_LEVELS:
            if rate == 0.0:
                simulation = clean
            else:
                simulation = simulate_expression(
                    meta,
                    seed,
                    spatial_basis="linear",
                    technical_noise={
                        "capture_efficiency": 1.0,
                        "ambient_rate": rate,
                        "seed": 40000 + seed,
                    },
                )
            realized_rate = float(
                simulation.omics_ambient_counts.to_numpy(dtype=float).sum()
                / len(simulation.meta)
            )
            calibration_rows.append(
                {
                    "seed": seed,
                    "component": "ambient",
                    "level": f"rate_{rate:g}",
                    "configured_value": rate,
                    "realized_value": realized_rate,
                }
            )
            if rate == 100.0:
                condition_simulations["Ambient 100"] = simulation

        for level, dropout_config in DROPOUT_LEVELS:
            if dropout_config is None:
                simulation = clean
            else:
                configured_dropout = dict(dropout_config)
                configured_dropout["seed"] = 50000 + seed
                simulation = simulate_expression(
                    meta,
                    seed,
                    spatial_basis="linear",
                    dropout=configured_dropout,
                )
            expected_probability = float(
                simulation.omics_dropout_probability.to_numpy(dtype=float).mean()
            )
            realized_probability = float(
                simulation.omics_dropout_mask.to_numpy(dtype=float).mean()
            )
            calibration_rows.append(
                {
                    "seed": seed,
                    "component": "dropout",
                    "level": level,
                    "configured_value": expected_probability,
                    "realized_value": realized_probability,
                }
            )
            if level == "High":
                condition_simulations["Dropout high"] = simulation

        combined_technical = dict(COMBINED_TECHNICAL_NOISE)
        combined_technical["seed"] = 60000 + seed
        combined_dropout = dict(COMBINED_DROPOUT)
        combined_dropout["seed"] = 70000 + seed
        combined = simulate_expression(
            meta,
            seed,
            spatial_basis="linear",
            technical_noise=combined_technical,
            dropout=combined_dropout,
        )
        condition_simulations["Combined"] = combined

        clean_median_library = float(
            np.median(clean.omics.to_numpy(dtype=float).sum(axis=1))
        )
        for condition in CONDITION_ORDER:
            metrics = expression_metrics(condition_simulations[condition], condition, seed)
            metrics["median_library_ratio"] = (
                float(metrics["median_library_size"]) / clean_median_library
            )
            observation_rows.append(metrics)

        if seed == REPRESENTATIVE_SEED:
            representative = representative_map(combined)
        print(
            f"Completed simulation seed {seed} ({SEEDS.index(seed) + 1}/{len(SEEDS)}).",
            flush=True,
        )

    if representative is None:
        raise RuntimeError("The representative map was not generated.")

    frames = {
        "effect": pd.DataFrame(effect_rows),
        "calibration": pd.DataFrame(calibration_rows),
        "observation": pd.DataFrame(observation_rows),
        "compatibility": pd.DataFrame(compatibility_rows),
        "truth": pd.DataFrame(truth_rows),
        "representative": representative,
    }
    checks = acceptance_checks(
        frames["effect"], frames["calibration"], frames["compatibility"]
    )
    frames["checks"] = checks
    if not checks["passed"].all():
        failures = checks.loc[~checks["passed"], ["check", "value", "criterion"]]
        raise RuntimeError(f"Predeclared acceptance checks failed:\n{failures.to_string(index=False)}")

    effect_seed_summary = (
        frames["effect"]
        .groupby(["seed", "basis", "configured_strength"], as_index=False)[
            ["latent_strength", "observed_strength"]
        ]
        .median()
    )
    summary_rows: list[dict[str, object]] = []
    summary_rows.extend(
        summarize_values(
            effect_seed_summary,
            ["basis", "configured_strength"],
            ["latent_strength", "observed_strength"],
            "direct_spatial_effect",
        )
    )
    summary_rows.extend(
        summarize_values(
            frames["calibration"],
            ["component", "level"],
            ["configured_value", "realized_value"],
            "observation_control_calibration",
        )
    )
    summary_rows.extend(
        summarize_values(
            frames["observation"],
            ["condition"],
            [
                "median_library_ratio",
                "zero_fraction",
                "mean_variance_spearman",
                "latent_recovery_pearson",
                "observed_spatial_slope",
            ],
            "observation_consequences",
        )
    )
    summary_rows.extend(
        {
            "section": "compatibility",
            "group": "all seeds",
            "metric": metric,
            "n_seeds": len(SEEDS),
            "median": float(frames["compatibility"][metric].mean()),
            "minimum": float(frames["compatibility"][metric].min()),
            "maximum": float(frames["compatibility"][metric].max()),
        }
        for metric in (
            "default_equals_disabled",
            "default_equals_zero_spatial",
            "default_disabled_mismatched_entries",
            "default_zero_spatial_mismatched_entries",
        )
    )
    frames["summary"] = pd.DataFrame(summary_rows)
    write_data(frames)
    write_configuration()
    return frames


def write_data(frames: dict[str, pd.DataFrame]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    frames["effect"].to_csv(DATA_DIR / "spatial_effect_recovery.tsv", sep="\t", index=False)
    frames["calibration"].to_csv(
        DATA_DIR / "observation_control_calibration.tsv", sep="\t", index=False
    )
    frames["observation"].to_csv(
        DATA_DIR / "observation_metrics.tsv", sep="\t", index=False
    )
    frames["compatibility"].to_csv(
        DATA_DIR / "backward_compatibility.tsv", sep="\t", index=False
    )
    frames["truth"].to_csv(DATA_DIR / "gene_truth.tsv", sep="\t", index=False)
    frames["representative"].to_csv(
        DATA_DIR / "representative_map.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    frames["checks"].to_csv(DATA_DIR / "acceptance_checks.tsv", sep="\t", index=False)
    frames["summary"].to_csv(DATA_DIR / "summary_metrics.tsv", sep="\t", index=False)


def write_configuration() -> None:
    configuration = {
        "analysis": "Supplementary Figure 12 optional spatial-expression and observation layers",
        "reviewer_comments": ["R2-M1", "R2-m5"],
        "related_comment": "R2-M3",
        "simspace_expected_version": SIMSPACE_EXPECTED_VERSION,
        "simspace_installed_version": package_version("simspace"),
        "simspace_source_commit": SIMSPACE_SOURCE_COMMIT,
        "adjacent_simspace_commit_at_run": adjacent_simspace_commit(),
        "adjacent_simspace_source_dirty_at_run": adjacent_simspace_source_dirty(),
        "python": sys.version,
        "platform": platform.platform(),
        "software": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
        },
        "seeds": list(SEEDS),
        "representative_seed": REPRESENTATIVE_SEED,
        "representative_gene_id": REPRESENTATIVE_GENE_ID,
        "representative_gene": REPRESENTATIVE_GENE,
        "representative_gene_selection": REPRESENTATIVE_GENE_SELECTION,
        "representative_gene_effect": {
            "basis": "linear",
            "configured_strength": 1.0,
            "direction": "Anti-diagonal",
        },
        "shape": list(SHAPE),
        "n_states": N_STATES,
        "n_genes": N_GENES,
        "mrf_sweeps": MRF_SWEEPS,
        "mrf_phi": MRF_PHI,
        "density_retention": DENSITY_RETENTION,
        "coordinate_jitter": COORDINATE_JITTER,
        "spatial_basis_families": list(SPATIAL_BASIS_ORDER),
        "spatial_basis_settings": {
            "linear_directions": {
                label: direction.tolist() for label, direction in LINEAR_DIRECTIONS.items()
            },
            "radial_center": [0.5, 0.5],
            "hotspot_center": [0.30, 0.70],
            "hotspot_bandwidth": 0.18,
            "structure_points": direct_spatial_config("structure_distance")[
                "structure_points"
            ],
            "coordinate_scaling": "unit_range",
            "center_basis": True,
        },
        "representative_spatial_direction": SPATIAL_DIRECTION.tolist(),
        "spatial_strengths": list(SPATIAL_STRENGTHS),
        "n_effect_genes_per_basis": N_EFFECT_GENES,
        "n_null_validation_genes_per_basis": N_NULL_VALIDATION_GENES,
        "capture_levels": list(CAPTURE_LEVELS),
        "ambient_rates_counts_per_cell": list(AMBIENT_LEVELS),
        "dropout_levels": [
            {"label": label, "configuration": configuration}
            for label, configuration in DROPOUT_LEVELS
        ],
        "combined_technical_noise": COMBINED_TECHNICAL_NOISE,
        "combined_dropout": COMBINED_DROPOUT,
        "generative_order": [
            "cell-type-conditioned Gamma means",
            "optional direct spatial log-mean effect",
            "Poisson biological count realization",
            "optional binomial capture thinning",
            "optional Poisson ambient background",
            "optional excess dropout",
        ],
        "statistical_unit": "simulation seed",
    }
    (DATA_DIR / "analysis_config.json").write_text(
        json.dumps(configuration, indent=2) + "\n", encoding="utf-8"
    )


def load_data() -> dict[str, pd.DataFrame]:
    required = {
        "effect": DATA_DIR / "spatial_effect_recovery.tsv",
        "calibration": DATA_DIR / "observation_control_calibration.tsv",
        "observation": DATA_DIR / "observation_metrics.tsv",
        "compatibility": DATA_DIR / "backward_compatibility.tsv",
        "truth": DATA_DIR / "gene_truth.tsv",
        "representative": DATA_DIR / "representative_map.tsv.gz",
        "checks": DATA_DIR / "acceptance_checks.tsv",
        "summary": DATA_DIR / "summary_metrics.tsv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing archived Supplementary Figure 12 data: " + ", ".join(missing))
    frames: dict[str, pd.DataFrame] = {}
    for name, path in required.items():
        frames[name] = pd.read_csv(path, sep="\t")
    return frames


def add_panel_label(fig: plt.Figure, spec, label: str) -> None:
    bounds = spec.get_position(fig)
    fig.text(bounds.x0 - 0.025, bounds.y1 + 0.012, label, fontsize=15, fontweight="bold")


def plot_seed_summary(
    axis: plt.Axes,
    frame: pd.DataFrame,
    x_column: str,
    y_column: str,
    color: str,
) -> None:
    rng = np.random.default_rng(912)
    x_values = np.sort(frame[x_column].unique().astype(float))
    for x_value in x_values:
        values = frame.loc[frame[x_column] == x_value, y_column].to_numpy(dtype=float)
        jitter = rng.normal(0.0, 0.012 * max(np.ptp(x_values), 1.0), size=len(values))
        axis.scatter(
            np.full(len(values), x_value) + jitter,
            values,
            s=16,
            color=color,
            alpha=0.35,
            edgecolors="none",
            rasterized=True,
        )
        median = float(np.median(values))
        axis.errorbar(
            x_value,
            median,
            yerr=[[median - values.min()], [values.max() - median]],
            fmt="o",
            color=color,
            markeredgecolor="white",
            markeredgewidth=0.7,
            markersize=5.5,
            linewidth=1.2,
            capsize=2,
            zorder=5,
        )


def draw_panel_a(fig: plt.Figure, spec, representative: pd.DataFrame) -> None:
    # Reserve identical map cells in the first row and place every legend or
    # colorbar in a separate second row. This prevents colorbars from shrinking
    # only the continuous maps and keeps all three spatial squares aligned.
    grid = spec.subgridspec(
        2,
        3,
        height_ratios=(1.0, 0.10),
        hspace=0.10,
        wspace=0.12,
    )
    axes = [fig.add_subplot(grid[0, index]) for index in range(3)]
    support_axes = [fig.add_subplot(grid[1, index]) for index in range(3)]
    state = representative["state"].to_numpy(dtype=int)
    axes[0].scatter(
        representative["col"],
        representative["row"],
        c=[STATE_COLORS[value] for value in state],
        s=4,
        linewidths=0,
        rasterized=True,
    )
    axes[0].set_title("Phenotype labels", fontsize=9)
    legend = [
        Line2D([0], [0], marker="o", linestyle="", color=STATE_COLORS[index], label=f"Type {index}", markersize=4)
        for index in range(N_STATES)
    ]
    support_axes[0].axis("off")
    support_axes[0].legend(
        handles=legend,
        frameon=False,
        fontsize=6.5,
        ncol=2,
        loc="center",
        columnspacing=0.6,
        handletextpad=0.2,
    )

    truth_limit = float(np.max(np.abs(representative["direct_log_fold_truth"])))
    truth_plot = axes[1].scatter(
        representative["col"],
        representative["row"],
        c=representative["direct_log_fold_truth"],
        cmap="coolwarm",
        vmin=-truth_limit,
        vmax=truth_limit,
        s=4,
        linewidths=0,
        rasterized=True,
    )
    axes[1].set_title("Direct log-fold truth", fontsize=9)
    colorbar = fig.colorbar(truth_plot, cax=support_axes[1], orientation="horizontal")
    colorbar.ax.tick_params(labelsize=6, length=2)

    observed_plot = axes[2].scatter(
        representative["col"],
        representative["row"],
        c=np.log1p(representative["observed_count"]),
        cmap="viridis",
        s=4,
        linewidths=0,
        rasterized=True,
    )
    axes[2].set_title(f"Observed {REPRESENTATIVE_GENE}", fontsize=9)
    colorbar = fig.colorbar(observed_plot, cax=support_axes[2], orientation="horizontal")
    colorbar.set_label("log(1 + count), combined artifacts", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6, length=2)

    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    add_panel_label(fig, spec, "A")


def draw_panel_b(fig: plt.Figure, spec, effect: pd.DataFrame) -> None:
    axis = fig.add_subplot(spec)
    seed_summary = (
        effect.groupby(["seed", "basis", "configured_strength"], as_index=False)[
            ["latent_strength", "observed_strength"]
        ]
        .median()
    )
    offsets = dict(zip(SPATIAL_BASIS_ORDER, (-0.036, -0.012, 0.012, 0.036), strict=True))
    rng = np.random.default_rng(912)
    for basis_name in SPATIAL_BASIS_ORDER:
        basis_data = seed_summary[seed_summary["basis"] == basis_name]
        color = SPATIAL_BASIS_COLORS[basis_name]
        marker = SPATIAL_BASIS_MARKERS[basis_name]
        for strength in SPATIAL_STRENGTHS:
            values = basis_data.loc[
                basis_data["configured_strength"] == strength, "observed_strength"
            ].to_numpy(dtype=float)
            x_position = float(strength + offsets[basis_name])
            jitter = rng.normal(0.0, 0.004, size=len(values))
            axis.scatter(
                np.full(len(values), x_position) + jitter,
                values,
                s=10,
                marker=marker,
                color=color,
                alpha=0.24,
                edgecolors="none",
                rasterized=True,
            )
            median = float(np.median(values))
            axis.errorbar(
                x_position,
                median,
                yerr=[[median - values.min()], [values.max() - median]],
                fmt=marker,
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.55,
                markersize=4.8,
                linewidth=0.9,
                capsize=1.6,
                zorder=5,
            )
    limits = (-1.15, 1.15)
    axis.plot(limits, limits, linestyle="--", color="#555555", linewidth=1.0)
    axis.scatter(
        SPATIAL_STRENGTHS,
        SPATIAL_STRENGTHS,
        marker="D",
        s=17,
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.8,
        label="Latent truth",
        zorder=6,
    )
    basis_handles = [
        Line2D(
            [0],
            [0],
            marker=SPATIAL_BASIS_MARKERS[basis_name],
            linestyle="",
            markerfacecolor=SPATIAL_BASIS_COLORS[basis_name],
            markeredgecolor="none",
            color="none",
            label=SPATIAL_BASIS_LABELS[basis_name],
            markersize=5,
        )
        for basis_name in SPATIAL_BASIS_ORDER
    ]
    basis_handles.extend(
        [
            Line2D(
                [0],
                [0],
                marker="D",
                linestyle="",
                markerfacecolor="none",
                markeredgecolor="#111111",
                color="none",
                label="Latent truth",
                markersize=4.5,
            ),
            Line2D(
                [0],
                [0],
                linestyle="--",
                color="#555555",
                label="Identity",
                linewidth=1.0,
            ),
        ]
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_xlabel("Configured log-mean coefficient")
    axis.set_ylabel("Recovered coefficient from observed counts")
    axis.set_title("Recovery across spatial bases", fontsize=11, loc="left")
    axis.legend(
        handles=basis_handles,
        frameon=False,
        fontsize=6.8,
        loc="upper left",
        ncol=2,
        columnspacing=0.7,
        handletextpad=0.35,
    )
    axis.grid(color="#DDDDDD", linewidth=0.5, alpha=0.7)
    add_panel_label(fig, spec, "B")


def draw_panel_c(fig: plt.Figure, spec, calibration: pd.DataFrame) -> None:
    grid = spec.subgridspec(2, 4, wspace=0.82, hspace=0.72)
    settings = (
        (
            "capture",
            "Capture efficiency",
            "Retained-count fraction",
            "#59A14F",
            grid[0, 0:2],
        ),
        (
            "ambient",
            "Ambient rate",
            "Added counts per cell",
            "#F28E2B",
            grid[0, 2:4],
        ),
        (
            "dropout",
            "Mean dropout probability",
            "Realized dropout fraction",
            "#B279A2",
            grid[1, 1:3],
        ),
    )
    for component, x_label, y_label, color, axis_spec in settings:
        axis = fig.add_subplot(axis_spec)
        component_data = calibration[calibration["component"] == component]
        plot_seed_summary(axis, component_data, "configured_value", "realized_value", color)
        values = np.concatenate(
            (
                component_data["configured_value"].to_numpy(dtype=float),
                component_data["realized_value"].to_numpy(dtype=float),
            )
        )
        lower = min(0.0, float(values.min()))
        upper = float(values.max())
        padding = max((upper - lower) * 0.08, 0.01)
        axis.plot(
            (lower - padding, upper + padding),
            (lower - padding, upper + padding),
            linestyle="--",
            color="#666666",
            linewidth=0.9,
        )
        axis.set_xlim(lower - padding, upper + padding)
        axis.set_ylim(lower - padding, upper + padding)
        axis.set_xlabel(x_label, fontsize=8)
        axis.set_ylabel(y_label, fontsize=8)
        axis.tick_params(labelsize=7)
        axis.grid(color="#E3E3E3", linewidth=0.45)
    bounds = spec.get_position(fig)
    fig.text(bounds.x0, bounds.y1 + 0.012, "Observation-control calibration", fontsize=11)
    add_panel_label(fig, spec, "C")


def condition_boxplot(
    axis: plt.Axes,
    observation: pd.DataFrame,
    metric: str,
    title: str,
    reference: float | None = None,
) -> None:
    values = [
        observation.loc[observation["condition"] == condition, metric].to_numpy(dtype=float)
        for condition in CONDITION_ORDER
    ]
    boxes = axis.boxplot(
        values,
        positions=np.arange(len(CONDITION_ORDER)),
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.1},
        whiskerprops={"color": "#555555", "linewidth": 0.8},
        capprops={"color": "#555555", "linewidth": 0.8},
    )
    for box, condition in zip(boxes["boxes"], CONDITION_ORDER, strict=True):
        box.set_facecolor(CONDITION_COLORS[condition])
        box.set_alpha(0.68)
        box.set_edgecolor("#555555")
        box.set_linewidth(0.7)
    rng = np.random.default_rng(203)
    for position, (condition, condition_values) in enumerate(zip(CONDITION_ORDER, values, strict=True)):
        jitter = rng.normal(0.0, 0.045, size=len(condition_values))
        axis.scatter(
            np.full(len(condition_values), position) + jitter,
            condition_values,
            s=8,
            color=CONDITION_COLORS[condition],
            alpha=0.32,
            edgecolors="none",
            rasterized=True,
        )
    if reference is not None:
        axis.axhline(reference, color="#777777", linestyle="--", linewidth=0.8)
    axis.set_xticks(np.arange(len(CONDITION_ORDER)))
    axis.set_xticklabels(("Clean", "Capture", "Ambient", "Dropout", "Combined"), rotation=28, ha="right")
    axis.set_title(title, fontsize=8.5)
    axis.tick_params(labelsize=7)
    axis.grid(axis="y", color="#E3E3E3", linewidth=0.45)


def draw_panel_d(fig: plt.Figure, spec, observation: pd.DataFrame) -> None:
    grid = spec.subgridspec(2, 2, wspace=0.35, hspace=0.62)
    settings = (
        ("median_library_ratio", "Median library-size ratio", 1.0),
        ("zero_fraction", "Zero fraction", None),
        ("mean_variance_spearman", "Mean–variance Spearman ρ", None),
        ("latent_recovery_pearson", "Latent–observed Pearson r", None),
    )
    for index, (metric, title, reference) in enumerate(settings):
        axis = fig.add_subplot(grid[index // 2, index % 2])
        condition_boxplot(axis, observation, metric, title, reference)
    bounds = spec.get_position(fig)
    title_y = bounds.y1 + 0.035
    fig.text(bounds.x0, title_y, "Consequences for observed count properties", fontsize=11)
    fig.text(bounds.x0 - 0.025, title_y, "D", fontsize=15, fontweight="bold")


def render_figure(frames: dict[str, pd.DataFrame]) -> None:
    sns.set_theme(style="ticks", context="paper")
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.8,
        }
    )
    figure = plt.figure(figsize=(13.0, 9.4), constrained_layout=False)
    outer = figure.add_gridspec(
        2,
        2,
        width_ratios=(1.25, 1.0),
        height_ratios=(0.85, 1.15),
        left=0.06,
        right=0.985,
        bottom=0.06,
        top=0.96,
        wspace=0.24,
        hspace=0.32,
    )
    draw_panel_a(figure, outer[0, 0], frames["representative"])
    draw_panel_b(figure, outer[0, 1], frames["effect"])
    draw_panel_c(figure, outer[1, 0], frames["calibration"])
    draw_panel_d(figure, outer[1, 1], frames["observation"])
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    png_path = SCRIPT_DIR / "SFig12.png"
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    shutil.copy2(png_path, EXAMPLE_DIR / png_path.name)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--render-only",
        # action="store_true",
        default=True,
        help="Render from archived TSV files without rerunning the simulations.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    frames = load_data() if arguments.render_only else generate_data()
    render_figure(frames)
    print(f"Wrote {SCRIPT_DIR / 'SFig12.png'}")
    if not arguments.render_only:
        print(f"All {len(frames['checks'])} predeclared acceptance checks passed.")


if __name__ == "__main__":
    main()
