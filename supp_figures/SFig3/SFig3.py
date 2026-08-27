"""Supplementary Figure 3: fixed spatial effects and observation layers.

This reproducible analysis addresses R2-M1 and R2-m5. It validates the
reference-free finite-basis fixed-effect model with a base/reference spatial
surface and treatment-coded cell-type deviations, then separately validates capture thinning,
ambient background, and mean-dependent dropout. The resulting conditional gene
truth could support a future separately prespecified known-truth SVG benchmark;
the existing R2-M3 concordance/stability analysis remains separate.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, distributions, version
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
from scipy.stats import chi2, pearsonr, spearmanr

import simspace as ss


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "Panel_A_D_data"
EXAMPLE_DIR = REPO_ROOT / "example_output" / "SFig3"

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
REFERENCE_STATE = 0
NON_REFERENCE_STATES = (1, 2, 3)
EFFECT_ROLES = ("overall", "state_1", "state_2", "state_3")
EFFECT_ROLE_LABELS = {
    "overall": "Base/reference f0",
    "state_1": "Type 1 deviation",
    "state_2": "Type 2 deviation",
    "state_3": "Type 3 deviation",
}
EFFECT_SCOPE_LABELS = {
    "overall": "Base/reference f0",
    "cell_type_specific": "Cell-type-specific fk",
}
N_EFFECT_REPLICATES = 3
N_NULLS_PER_ROLE = 6
DETECTION_ALPHA = 0.05
SPATIAL_DIRECTION = np.asarray((1.0, 1.0), dtype=float) / np.sqrt(2.0)
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
SIMSPACE_EXPECTED_VERSION = "0.4.1"
REPRESENTATIVE_SEED = 0
REPRESENTATIVE_GENE_SELECTION = (
    "For the illustrative seed-0 map only, rank the nine linear genes with a "
    "cell-type-specific deviation of +1 using a deterministic, count-blind "
    "visualization score. The score combines configured truth span, expected "
    "post-artifact visibility, Poisson working contrast and design-information "
    "proxies, direct-effect attribution, and target-to-background prominence. "
    "Break score ties by GeneID. Realized counts do not enter selection, and the "
    "display choice does not enter the coefficient-recovery or likelihood-ratio "
    "summaries."
)
REPRESENTATIVE_SCORE_WEIGHTS = {
    "truth_q90_span": 0.22,
    "poisson_working_contrast_proxy": 0.25,
    "direct_effect_fraction_of_log_count_sd": 0.23,
    "expected_visible_fraction": 0.12,
    "sqrt_poisson_working_information_proxy": 0.10,
    "target_high_to_nontarget_q95": 0.08,
}

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


def visible_package_distributions(package: str) -> list[dict[str, str]]:
    """Record every matching distribution visible on the effective sys.path."""
    return [
        {
            "version": distribution.version,
            "metadata_path": str(distribution._path),
        }
        for distribution in distributions(name=package)
    ]


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    """Declare balanced base/reference, deviation, and null-control genes."""
    if basis_name not in SPATIAL_BASIS_ORDER:
        raise ValueError(f"Unknown spatial basis: {basis_name}")

    rows: list[dict[str, object]] = []
    gene_id = 0
    for role_index, role in enumerate(EFFECT_ROLES):
        target_state = None if role == "overall" else int(role.removeprefix("state_"))
        effect_scope = "overall" if role == "overall" else "cell_type_specific"
        for strength_index, strength in enumerate(NONZERO_SPATIAL_STRENGTHS):
            for replicate in range(N_EFFECT_REPLICATES):
                if basis_name == "linear":
                    direction_name = "Diagonal"
                    direction = SPATIAL_DIRECTION
                else:
                    direction_name = "Not applicable"
                    direction = np.asarray((np.nan, np.nan), dtype=float)
                rows.append(
                    {
                        "GeneID": gene_id,
                        "Gene": f"Gene_{gene_id}",
                        "configured_strength": float(strength),
                        "effect_role": role,
                        "effect_role_label": EFFECT_ROLE_LABELS[role],
                        "effect_scope": effect_scope,
                        "effect_scope_label": EFFECT_SCOPE_LABELS[effect_scope],
                        "target_state": target_state,
                        "direction": direction_name,
                        "direction_row": float(direction[0]),
                        "direction_col": float(direction[1]),
                        "replicate": replicate,
                        "is_direct": True,
                        "is_null": False,
                    }
                )
                gene_id += 1

    if len(rows) != N_EFFECT_GENES:
        raise RuntimeError("The declared spatial-effect design has the wrong size.")

    for role_index, role in enumerate(EFFECT_ROLES):
        target_state = None if role == "overall" else int(role.removeprefix("state_"))
        effect_scope = "overall" if role == "overall" else "cell_type_specific"
        for replicate in range(N_NULLS_PER_ROLE):
            if basis_name == "linear":
                direction_name = "Diagonal"
                direction = SPATIAL_DIRECTION
            else:
                direction_name = "Not applicable"
                direction = np.asarray((np.nan, np.nan), dtype=float)
            rows.append(
                {
                    "GeneID": gene_id,
                    "Gene": f"Gene_{gene_id}",
                    "configured_strength": 0.0,
                    "effect_role": role,
                    "effect_role_label": EFFECT_ROLE_LABELS[role],
                    "effect_scope": effect_scope,
                    "effect_scope_label": EFFECT_SCOPE_LABELS[effect_scope],
                    "target_state": target_state,
                    "direction": direction_name,
                    "direction_row": float(direction[0]),
                    "direction_col": float(direction[1]),
                    "replicate": replicate,
                    "is_direct": False,
                    "is_null": True,
                }
            )
            gene_id += 1

    design = pd.DataFrame(rows)
    nonzero_counts = design.loc[design["is_direct"]].groupby(
        ["effect_role", "configured_strength"]
    ).size()
    null_counts = design.loc[design["is_null"]].groupby("effect_role").size()
    if not (nonzero_counts == N_EFFECT_REPLICATES).all():
        raise RuntimeError("Nonzero effects are not balanced across strengths and roles.")
    if not (null_counts == N_NULLS_PER_ROLE).all():
        raise RuntimeError("Null controls are not balanced across test roles.")
    if len(design) != N_EFFECT_GENES + N_NULL_VALIDATION_GENES:
        raise RuntimeError("The declared validation design has the wrong total size.")
    return design


def configured_coefficient_vector(effect_row: object, basis_name: str) -> np.ndarray:
    """Return the vector attached to the configured base or deviation block."""
    strength = float(effect_row.configured_strength)
    if basis_name == "linear":
        return strength * np.asarray(
            [float(effect_row.direction_row), float(effect_row.direction_col)],
            dtype=float,
        )
    return np.asarray([strength], dtype=float)


def direct_spatial_config(
    basis_name: str,
    *,
    zero_coefficients: bool = False,
) -> dict[str, object]:
    """Build the finite-basis treatment-coded configuration for one basis."""
    design = validation_effect_design(basis_name)
    direct = design[design["is_direct"]].copy()
    n_basis = 2 if basis_name == "linear" else 1
    overall_coefficients: dict[str, list[float]] = {}
    cell_type_coefficients: dict[str, dict[int, list[float]]] = {}
    for row in direct.itertuples(index=False):
        coefficient = configured_coefficient_vector(row, basis_name)
        if zero_coefficients:
            coefficient = np.zeros(n_basis, dtype=float)
        if row.effect_role == "overall":
            overall_coefficients[row.Gene] = coefficient.tolist()
        else:
            overall_coefficients[row.Gene] = np.zeros(n_basis, dtype=float).tolist()
            cell_type_coefficients[row.Gene] = {
                int(row.target_state): coefficient.tolist()
            }

    config: dict[str, object] = {
        "genes": direct["Gene"].tolist(),
        "basis": basis_name,
        "overall_coefficients": overall_coefficients,
        "cell_type_coefficients": cell_type_coefficients,
        "reference_state": REFERENCE_STATE,
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


def apply_observation_variant(
    clean: ss.SimSpace,
    *,
    technical_noise: dict[str, object] | None = None,
    dropout: dict[str, object] | None = None,
) -> ss.SimSpace:
    """Apply package observation functions to one fixed biological draw."""
    simulation = copy.copy(clean)
    resolved_technical = technical_noise or {}
    observation_seed = int(resolved_technical.get("seed", clean.seed + 10000))
    capture_counts, capture_efficiency = ss.omics.applyCaptureNoise(
        clean.omics_baseline_counts,
        capture_efficiency=resolved_technical.get("capture_efficiency", 1.0),
        seed=observation_seed,
    )
    ambient_counts, ambient_additions, ambient_rate, ambient_profile = (
        ss.omics.applyAmbientNoise(
            capture_counts,
            ambient_rate=resolved_technical.get("ambient_rate", 0.0),
            ambient_profile=resolved_technical.get("ambient_profile"),
            seed=observation_seed + 1,
        )
    )
    dropout_seed = int((dropout or {}).get("seed", clean.seed + 20000))
    observed, dropout_mask, dropout_probability = ss.omics.applyDropout(
        ambient_counts,
        dropout=dropout,
        latent_mean=clean.omics_latent_mean,
        seed=dropout_seed,
    )
    simulation.omics = observed
    simulation.omics_capture_counts = capture_counts
    simulation.omics_ambient_counts = ambient_additions
    simulation.omics_dropout_mask = dropout_mask
    simulation.omics_dropout_probability = dropout_probability
    simulation.omics_observed_counts = observed.copy()
    simulation.omics_observation_truth = {
        "capture_efficiency": capture_efficiency,
        "ambient_rate": ambient_rate,
        "ambient_profile": ambient_profile,
        "observation_seed": observation_seed,
        "dropout_seed": dropout_seed,
    }
    return simulation


def spatial_basis_matrix(simulation: ss.SimSpace, basis_name: str) -> np.ndarray:
    """Extract the complete configured basis matrix from stored simulation truth."""
    basis_columns = {
        "linear": ["row_linear", "col_linear"],
        "radial": ["radial_distance"],
        "hotspot": ["hotspot_basis"],
        "structure_distance": ["structure_distance"],
    }[basis_name]
    return simulation.omics_spatial_design[basis_columns].to_numpy(dtype=float)


def full_fixed_effect_design(
    meta: pd.DataFrame,
    basis: np.ndarray,
) -> tuple[np.ndarray, dict[str, slice]]:
    """Build intercept, phenotype main effects, B, and phenotype-by-B terms."""
    states = meta["state"].to_numpy(dtype=int)
    observed_states = set(np.unique(states).tolist())
    expected_states = {REFERENCE_STATE, *NON_REFERENCE_STATES}
    if observed_states != expected_states:
        raise ValueError(
            f"Expected phenotype states {sorted(expected_states)}; found {sorted(observed_states)}."
        )
    main_effects = np.column_stack([states == state for state in NON_REFERENCE_STATES]).astype(float)
    n_basis = basis.shape[1]
    columns = [np.ones((len(meta), 1), dtype=float), main_effects, basis]
    blocks = {"overall": slice(1 + len(NON_REFERENCE_STATES), 1 + len(NON_REFERENCE_STATES) + n_basis)}
    start = blocks["overall"].stop
    for state_index, state in enumerate(NON_REFERENCE_STATES):
        columns.append(main_effects[:, [state_index]] * basis)
        blocks[f"state_{state}"] = slice(start, start + n_basis)
        start += n_basis
    return np.column_stack(columns), blocks


def poisson_log_likelihood(
    counts: np.ndarray,
    design: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    linear_predictor = np.clip(design @ coefficients, -30.0, 30.0)
    return np.sum(counts * linear_predictor - np.exp(linear_predictor), axis=0)


def poisson_log_link_fit_batch(
    counts: np.ndarray,
    design: np.ndarray,
    *,
    initial: np.ndarray | None = None,
    tolerance: float = 1e-8,
    max_iterations: int = 60,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit independent Poisson GLMs sharing one design by Fisher scoring."""
    counts = np.asarray(counts, dtype=float)
    if counts.ndim == 1:
        counts = counts[:, None]
    design = np.asarray(design, dtype=float)
    if counts.shape[0] != design.shape[0]:
        raise ValueError("Counts and fixed-effect design have different cell counts.")
    if np.any(counts < 0) or not np.isfinite(counts).all():
        raise ValueError("Poisson counts must be finite and non-negative.")

    n_coefficients = design.shape[1]
    n_genes = counts.shape[1]
    coefficients = (
        np.zeros((n_coefficients, n_genes), dtype=float)
        if initial is None
        else np.asarray(initial, dtype=float).copy()
    )
    if coefficients.shape != (n_coefficients, n_genes):
        raise ValueError("Initial coefficients do not match the design and gene count.")

    converged = np.zeros(n_genes, dtype=bool)
    iteration_count = np.full(n_genes, max_iterations, dtype=int)
    diagonal = np.arange(n_coefficients)
    for iteration in range(1, max_iterations + 1):
        linear_predictor = np.clip(design @ coefficients, -30.0, 30.0)
        fitted_mean = np.exp(linear_predictor)
        score = design.T @ (counts - fitted_mean)
        information = np.einsum(
            "ni,nm,nj->mij",
            design,
            fitted_mean,
            design,
            optimize=True,
        )
        diagonal_scale = np.maximum(
            information[:, diagonal, diagonal].mean(axis=1),
            1.0,
        )
        information[:, diagonal, diagonal] += diagonal_scale[:, None] * 1e-10
        try:
            step = np.linalg.solve(information, score.T[:, :, None])[:, :, 0].T
        except np.linalg.LinAlgError:
            step = np.column_stack(
                [
                    np.linalg.lstsq(information[index], score[:, index], rcond=None)[0]
                    for index in range(n_genes)
                ]
            )
        step[:, converged] = 0.0
        step = np.clip(step, -5.0, 5.0)

        current_likelihood = np.sum(
            counts * linear_predictor - fitted_mean,
            axis=0,
        )
        scales = np.ones(n_genes, dtype=float)
        candidate = coefficients + step
        candidate_likelihood = poisson_log_likelihood(counts, design, candidate)
        for _ in range(14):
            worse = candidate_likelihood < current_likelihood - 1e-8
            if not np.any(worse):
                break
            scales[worse] *= 0.5
            candidate[:, worse] = (
                coefficients[:, worse] + step[:, worse] * scales[worse][None, :]
            )
            candidate_likelihood[worse] = poisson_log_likelihood(
                counts[:, worse],
                design,
                candidate[:, worse],
            )
        worse = candidate_likelihood < current_likelihood - 1e-8
        if np.any(worse):
            failed = np.flatnonzero(worse).tolist()
            raise RuntimeError(
                "Poisson Fisher-scoring line search failed for batch column(s) "
                f"{failed}."
            )
        delta = np.max(np.abs(candidate - coefficients), axis=0)
        coefficients = candidate
        newly_converged = (~converged) & (delta < tolerance)
        iteration_count[newly_converged] = iteration
        converged |= newly_converged
        if converged.all():
            break

    likelihood = poisson_log_likelihood(counts, design, coefficients)
    raw_linear_predictor = design @ coefficients
    if np.max(np.abs(raw_linear_predictor)) >= 29.0:
        raise RuntimeError(
            "A fitted Poisson linear predictor approached the numerical clipping "
            "boundary; rescale the design before interpreting this fit."
        )
    return coefficients, likelihood, converged, iteration_count


def full_model_initial_coefficients(
    counts: np.ndarray,
    meta: pd.DataFrame,
    n_coefficients: int,
) -> np.ndarray:
    """Initialize the full GLM at phenotype-specific empirical means."""
    states = meta["state"].to_numpy(dtype=int)
    means = {
        state: np.clip(counts[states == state].mean(axis=0), 1e-8, None)
        for state in (REFERENCE_STATE, *NON_REFERENCE_STATES)
    }
    initial = np.zeros((n_coefficients, counts.shape[1]), dtype=float)
    initial[0] = np.log(means[REFERENCE_STATE])
    for column, state in enumerate(NON_REFERENCE_STATES, start=1):
        initial[column] = np.log(means[state]) - np.log(means[REFERENCE_STATE])
    return initial


def fit_spatial_effects(
    simulation: ss.SimSpace,
    basis_name: str,
) -> pd.DataFrame:
    """Fit the full finite-basis shared-plus-deviation design and block LRTs."""
    effect_design = validation_effect_design(basis_name).reset_index(drop=True)
    basis = spatial_basis_matrix(simulation, basis_name)
    fixed_design, blocks = full_fixed_effect_design(simulation.meta, basis)
    design_rank = int(np.linalg.matrix_rank(fixed_design))
    design_condition_number = float(np.linalg.cond(fixed_design))
    gene_names = effect_design["Gene"].tolist()
    counts = simulation.omics[gene_names].to_numpy(dtype=float)
    latent_mean = simulation.omics_latent_mean[gene_names].to_numpy(dtype=float)
    initial = full_model_initial_coefficients(counts, simulation.meta, fixed_design.shape[1])
    estimated, full_likelihood, full_converged, full_iterations = (
        poisson_log_link_fit_batch(counts, fixed_design, initial=initial)
    )
    latent_coefficients = np.linalg.lstsq(
        fixed_design,
        np.log(latent_mean),
        rcond=None,
    )[0]

    reduced_likelihood = np.full(len(effect_design), np.nan, dtype=float)
    reduced_converged = np.zeros(len(effect_design), dtype=bool)
    reduced_iterations = np.zeros(len(effect_design), dtype=int)
    for role in EFFECT_ROLES:
        positions = np.flatnonzero(effect_design["effect_role"].to_numpy() == role)
        dropped = np.arange(blocks[role].start, blocks[role].stop)
        retained = np.setdiff1d(np.arange(fixed_design.shape[1]), dropped)
        reduced_estimated, role_likelihood, role_converged, role_iterations = (
            poisson_log_link_fit_batch(
                counts[:, positions],
                fixed_design[:, retained],
                initial=estimated[retained][:, positions],
            )
        )
        del reduced_estimated
        reduced_likelihood[positions] = role_likelihood
        reduced_converged[positions] = role_converged
        reduced_iterations[positions] = role_iterations

    likelihood_difference = full_likelihood - reduced_likelihood
    likelihood_ratio = np.maximum(2.0 * likelihood_difference, 0.0)
    n_basis = basis.shape[1]
    p_values = chi2.sf(likelihood_ratio, df=n_basis)
    rows: list[dict[str, object]] = []
    for position, effect in enumerate(effect_design.itertuples(index=False)):
        block = blocks[str(effect.effect_role)]
        configured_vector = configured_coefficient_vector(effect, basis_name)
        latent_vector = latent_coefficients[block, position]
        estimated_vector = estimated[block, position]
        spatial_start = 1 + len(NON_REFERENCE_STATES)
        expected_spatial_vector = np.zeros(
            fixed_design.shape[1] - spatial_start,
            dtype=float,
        )
        target_slice = slice(block.start - spatial_start, block.stop - spatial_start)
        expected_spatial_vector[target_slice] = configured_vector
        latent_spatial_vector_error = float(
            np.max(
                np.abs(
                    latent_coefficients[spatial_start:, position]
                    - expected_spatial_vector
                )
            )
        )
        projection_direction = (
            np.asarray([float(effect.direction_row), float(effect.direction_col)])
            if basis_name == "linear"
            else np.ones(1, dtype=float)
        )
        latent_strength = float(latent_vector @ projection_direction)
        observed_strength = float(estimated_vector @ projection_direction)
        orthogonal_magnitude = float(
            np.linalg.norm(estimated_vector - observed_strength * projection_direction)
            if basis_name == "linear"
            else 0.0
        )
        rows.append(
            {
                "seed": int(simulation.seed - 10000),
                "basis": basis_name,
                "basis_label": SPATIAL_BASIS_LABELS[basis_name],
                "basis_dimension": n_basis,
                "GeneID": int(effect.GeneID),
                "Gene": str(effect.Gene),
                "is_direct": bool(effect.is_direct),
                "is_null": bool(effect.is_null),
                "effect_role": str(effect.effect_role),
                "effect_role_label": str(effect.effect_role_label),
                "effect_scope": str(effect.effect_scope),
                "effect_scope_label": str(effect.effect_scope_label),
                "reference_state": REFERENCE_STATE,
                "target_state": effect.target_state,
                "direction": str(effect.direction),
                "configured_strength": float(effect.configured_strength),
                "configured_coefficient_0": float(configured_vector[0]),
                "configured_coefficient_1": (
                    float(configured_vector[1]) if n_basis == 2 else np.nan
                ),
                "latent_coefficient_0": float(latent_vector[0]),
                "latent_coefficient_1": (
                    float(latent_vector[1]) if n_basis == 2 else np.nan
                ),
                "estimated_coefficient_0": float(estimated_vector[0]),
                "estimated_coefficient_1": (
                    float(estimated_vector[1]) if n_basis == 2 else np.nan
                ),
                "latent_strength": latent_strength,
                "observed_strength": observed_strength,
                "maximum_latent_spatial_vector_error": latent_spatial_vector_error,
                "orthogonal_estimate_magnitude": orthogonal_magnitude,
                "full_log_likelihood": float(full_likelihood[position]),
                "reduced_log_likelihood": float(reduced_likelihood[position]),
                "likelihood_ratio": float(likelihood_ratio[position]),
                "full_minus_reduced_log_likelihood": float(
                    likelihood_difference[position]
                ),
                "lrt_df": n_basis,
                "p_value": float(p_values[position]),
                "detection_alpha": DETECTION_ALPHA,
                "detected": bool(p_values[position] < DETECTION_ALPHA),
                "full_fit_converged": bool(full_converged[position]),
                "reduced_fit_converged": bool(reduced_converged[position]),
                "full_fit_iterations": int(full_iterations[position]),
                "reduced_fit_iterations": int(reduced_iterations[position]),
                "fixed_design_rank": design_rank,
                "fixed_design_columns": fixed_design.shape[1],
                "fixed_design_condition_number": design_condition_number,
                "n_cells": len(simulation.meta),
                "n_target_cells": int(
                    np.sum(
                        simulation.meta["state"].to_numpy(dtype=int)
                        == REFERENCE_STATE
                    )
                    if effect.effect_role == "overall"
                    else np.sum(
                        simulation.meta["state"].to_numpy(dtype=int)
                        == int(effect.target_state)
                    )
                ),
            }
        )
    return pd.DataFrame(rows)


def target_block_strength(
    simulation: ss.SimSpace,
    effect_row: object,
    basis_name: str = "linear",
) -> tuple[float, float]:
    """Estimate one configured block without an oracle phenotype-mean offset."""
    basis = spatial_basis_matrix(simulation, basis_name)
    fixed_design, blocks = full_fixed_effect_design(simulation.meta, basis)
    gene_name = str(effect_row.Gene)
    counts = simulation.omics[[gene_name]].to_numpy(dtype=float)
    initial = full_model_initial_coefficients(counts, simulation.meta, fixed_design.shape[1])
    estimated, _, converged, _ = poisson_log_link_fit_batch(
        counts,
        fixed_design,
        initial=initial,
    )
    if not converged[0]:
        raise RuntimeError(f"Full fixed-effect fit did not converge for {gene_name}.")
    latent = np.linalg.lstsq(
        fixed_design,
        np.log(simulation.omics_latent_mean[[gene_name]].to_numpy(dtype=float)),
        rcond=None,
    )[0]
    block = blocks[str(effect_row.effect_role)]
    direction = np.asarray(
        [float(effect_row.direction_row), float(effect_row.direction_col)],
        dtype=float,
    )
    return float(latent[block, 0] @ direction), float(estimated[block, 0] @ direction)


def expression_metrics(
    simulation: ss.SimSpace,
    condition: str,
    seed: int,
    representative_effect: object,
) -> dict[str, object]:
    latent_slope, observed_slope = target_block_strength(
        simulation,
        representative_effect,
    )
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
        is_validation_gene = gene_id in effect_design.index
        is_direct = is_validation_gene and bool(effect_design.loc[gene_id, "is_direct"])
        is_null = is_validation_gene and bool(effect_design.loc[gene_id, "is_null"])
        if is_validation_gene:
            effect = effect_design.loc[gene_id]
            tested_vector = configured_coefficient_vector(effect, "linear")
            coefficient_row = float(tested_vector[0])
            coefficient_col = float(tested_vector[1])
        else:
            effect = None
            coefficient_row = 0.0
            coefficient_col = 0.0
        if is_direct:
            log_effect = np.log(
                simulation.omics_latent_mean[f"Gene_{gene_id}"].to_numpy(dtype=float)
                / baseline[f"Gene_{gene_id}"].to_numpy(dtype=float)
            )
        else:
            log_effect = np.zeros(len(simulation.meta), dtype=float)
        effect_role = str(effect["effect_role"]) if effect is not None else "none"
        target_state = effect["target_state"] if effect is not None else np.nan
        overall_vector = (
            np.asarray([coefficient_row, coefficient_col], dtype=float)
            if is_direct and effect_role == "overall"
            else np.zeros(2, dtype=float)
        )
        interaction_vectors = {state: np.zeros(2, dtype=float) for state in NON_REFERENCE_STATES}
        if is_direct and effect_role != "overall":
            interaction_vectors[int(target_state)] = np.asarray(
                [coefficient_row, coefficient_col],
                dtype=float,
            )
        row: dict[str, object] = {
            "seed": seed,
            "GeneID": gene_id,
            "Gene": f"Gene_{gene_id}",
            "Marker": int(gene["Marker"]),
            "truth_scope": (
                "base_reference_spatial_fixed_effect"
                if is_direct and effect_role == "overall"
                else "cell_type_specific_spatial_deviation"
                if is_direct
                else "spatial_null_control"
                if is_null
                else "cell_type_only"
            ),
            "conditional_spatial_truth": is_direct,
            "is_null_control": is_null,
            "spatial_basis": "linear" if is_direct else "none",
            "effect_role": effect_role,
            "effect_role_label": (
                str(effect["effect_role_label"]) if effect is not None else "Not tested"
            ),
            "effect_scope": (
                str(effect["effect_scope"]) if effect is not None else "none"
            ),
            "reference_state": REFERENCE_STATE,
            "target_state": target_state,
            "configured_strength": (
                float(effect["configured_strength"]) if effect is not None else 0.0
            ),
            "direction": str(effect["direction"]) if effect is not None else "none",
            "coefficient_row": coefficient_row,
            "coefficient_col": coefficient_col,
            "overall_coefficient_0": float(overall_vector[0]),
            "overall_coefficient_1": float(overall_vector[1]),
            "cell_type_1_coefficient_0": float(interaction_vectors[1][0]),
            "cell_type_1_coefficient_1": float(interaction_vectors[1][1]),
            "cell_type_2_coefficient_0": float(interaction_vectors[2][0]),
            "cell_type_2_coefficient_1": float(interaction_vectors[2][1]),
            "cell_type_3_coefficient_0": float(interaction_vectors[3][0]),
            "cell_type_3_coefficient_1": float(interaction_vectors[3][1]),
            "min_log_effect": float(log_effect.min()),
            "max_log_effect": float(log_effect.max()),
        }
        for column in type_columns:
            row[column] = float(gene[column])
        rows.append(row)
    return rows


def select_representative_effect(
    clean: ss.SimSpace,
) -> tuple[pd.Series, pd.DataFrame]:
    """Rank display candidates using configured truth and expected counts only."""
    design = validation_effect_design("linear")
    candidates = design[
        (design["is_direct"])
        & (design["effect_scope"] == "cell_type_specific")
        & (design["configured_strength"] == 1.0)
    ]
    gene_meta = clean.gene_meta.set_index("GeneID")
    baseline = ss.omics.buildOmicsMean(clean.gene_meta, clean.meta)
    latent = clean.omics_latent_mean
    latent_values = latent.to_numpy(dtype=float)
    capture_efficiency = float(COMBINED_TECHNICAL_NOISE["capture_efficiency"])
    ambient_rate = float(COMBINED_TECHNICAL_NOISE["ambient_rate"])
    dropout_intercept = float(COMBINED_DROPOUT["intercept"])
    dropout_slope = float(COMBINED_DROPOUT["slope"])

    # The default ambient profile is proportional to captured feature totals + 1.
    # Replace the random capture draw with its expectation so selection remains
    # independent of the realized Poisson and observation-layer random numbers.
    expected_capture_totals = capture_efficiency * latent_values.sum(axis=0)
    ambient_profile = pd.Series(
        expected_capture_totals + 1.0,
        index=latent.columns,
        dtype=float,
    )
    ambient_profile /= ambient_profile.sum()

    basis_score = (
        clean.omics_spatial_design[["row_linear", "col_linear"]].to_numpy(dtype=float)
        @ SPATIAL_DIRECTION
    )
    basis_matrix = spatial_basis_matrix(clean, "linear")
    full_design, full_blocks = full_fixed_effect_design(clean.meta, basis_matrix)
    state = clean.meta["state"].to_numpy(dtype=int)
    scores: list[dict[str, float | int | str]] = []
    for effect in candidates.itertuples(index=False):
        gene_id = int(effect.GeneID)
        gene = str(effect.Gene)
        target_state = int(effect.target_state)
        target_mask = state == target_state
        latent_mean = latent[gene].to_numpy(dtype=float)
        baseline_mean = baseline[gene].to_numpy(dtype=float)
        truth = np.log(latent_mean / baseline_mean)

        ambient_mean = ambient_rate * float(ambient_profile.loc[gene])
        pre_dropout_mean = capture_efficiency * latent_mean + ambient_mean
        dropout_probability = 1.0 / (
            1.0
            + np.exp(
                -(
                    dropout_intercept
                    + dropout_slope * np.log1p(latent_mean)
                )
            )
        )
        expected_observed = (1.0 - dropout_probability) * pre_dropout_mean
        expected_visible_probability = (1.0 - dropout_probability) * (
            1.0 - np.exp(-pre_dropout_mean)
        )

        pre_dropout_baseline = capture_efficiency * baseline_mean + ambient_mean
        baseline_dropout_probability = 1.0 / (
            1.0
            + np.exp(
                -(
                    dropout_intercept
                    + dropout_slope * np.log1p(baseline_mean)
                )
            )
        )
        expected_observed_baseline = (
            1.0 - baseline_dropout_probability
        ) * pre_dropout_baseline

        target_truth = truth[target_mask]
        low_cut, high_cut = np.quantile(target_truth, [0.20, 0.80])
        low_mask = target_mask & (truth <= low_cut)
        high_mask = target_mask & (truth >= high_cut)
        low_mean = float(expected_observed[low_mask].mean())
        high_mean = float(expected_observed[high_mask].mean())
        expected_contrast = abs(high_mean - low_mean)
        poisson_working_contrast_proxy = expected_contrast / np.sqrt(
            high_mean + low_mean + 1e-12
        )

        # Poisson working design-information proxy for the configured diagonal
        # direction after projecting out all full-model nuisance columns,
        # including the target block's orthogonal direction. This is a display
        # score, not literal Fisher information for the dropout-contaminated law.
        target_block = full_blocks[str(effect.effect_role)]
        target_columns = np.arange(target_block.start, target_block.stop)
        orthogonal_direction = np.asarray(
            (-SPATIAL_DIRECTION[1], SPATIAL_DIRECTION[0]),
            dtype=float,
        )
        target_score = target_mask.astype(float) * basis_score
        orthogonal_score = target_mask.astype(float) * (
            basis_matrix @ orthogonal_direction
        )
        nuisance = np.column_stack(
            [
                np.delete(full_design, target_columns, axis=1),
                orthogonal_score,
            ]
        )
        sqrt_weights = np.sqrt(np.maximum(expected_observed, 1e-10))
        weighted_target_score = sqrt_weights * target_score
        weighted_nuisance = sqrt_weights[:, None] * nuisance
        nuisance_fit = np.linalg.lstsq(
            weighted_nuisance,
            weighted_target_score,
            rcond=None,
        )[0]
        efficient_residual = (
            weighted_target_score - weighted_nuisance @ nuisance_fit
        )

        direct_log_count_delta = np.log1p(expected_observed) - np.log1p(
            expected_observed_baseline
        )
        expected_log_count = np.log1p(expected_observed)
        scores.append(
            {
                "GeneID": gene_id,
                "Gene": gene,
                "basis": "linear",
                "effect_scope": str(effect.effect_scope),
                "effect_role": str(effect.effect_role),
                "effect_role_label": str(effect.effect_role_label),
                "configured_strength": float(effect.configured_strength),
                "direction": str(effect.direction),
                "target_state": target_state,
                "replicate": int(effect.replicate),
                "n_target_cells": int(target_mask.sum()),
                "target_baseline_mean": float(
                    gene_meta.loc[gene_id, f"Type_{target_state}"]
                ),
                "truth_q90_span": float(
                    np.quantile(target_truth, 0.95)
                    - np.quantile(target_truth, 0.05)
                ),
                "expected_visible_fraction": float(
                    expected_visible_probability[target_mask].mean()
                ),
                "expected_lower_quintile_count": low_mean,
                "expected_upper_quintile_count": high_mean,
                "expected_count_contrast": expected_contrast,
                "poisson_working_contrast_proxy": float(
                    poisson_working_contrast_proxy
                ),
                "sqrt_poisson_working_information_proxy": float(
                    np.linalg.norm(efficient_residual)
                ),
                "direct_effect_fraction_of_log_count_sd": float(
                    np.std(direct_log_count_delta)
                    / (np.std(expected_log_count) + 1e-12)
                ),
                "target_high_to_nontarget_q95": float(
                    max(high_mean, low_mean)
                    / (np.quantile(expected_observed[~target_mask], 0.95) + 1e-12)
                ),
            }
        )
    ranking = pd.DataFrame(scores)
    score_columns = list(REPRESENTATIVE_SCORE_WEIGHTS)
    if not np.isfinite(ranking[score_columns].to_numpy(dtype=float)).all():
        raise RuntimeError("Representative-gene display scores must all be finite.")
    for component in score_columns:
        ranking[f"{component}_percentile"] = ranking[component].rank(
            method="average",
            pct=True,
        )
    ranking["visualization_score"] = 100.0 * sum(
        weight * ranking[f"{component}_percentile"]
        for component, weight in REPRESENTATIVE_SCORE_WEIGHTS.items()
    )
    ranking = ranking.sort_values(
        ["visualization_score", "GeneID"],
        ascending=[False, True],
    ).reset_index(drop=True)
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    selected_gene_id = int(ranking.iloc[0]["GeneID"])
    selected = design.loc[design["GeneID"] == selected_gene_id].iloc[0].copy()
    for column, value in ranking.iloc[0].items():
        if column not in selected.index:
            selected[column] = value
    selected["selection_rank"] = int(ranking.iloc[0]["rank"])
    return selected, ranking


def representative_map(simulation: ss.SimSpace, effect: object) -> pd.DataFrame:
    selected_gene_id = int(effect.GeneID)
    representative_gene = str(effect.Gene)
    baseline = ss.omics.buildOmicsMean(simulation.gene_meta, simulation.meta)
    latent = simulation.omics_latent_mean[representative_gene].to_numpy(dtype=float)
    base = baseline[representative_gene].to_numpy(dtype=float)
    return pd.DataFrame(
        {
            "state": simulation.meta["state"].astype(int).to_numpy(),
            "row": simulation.meta["row"].to_numpy(dtype=float),
            "col": simulation.meta["col"].to_numpy(dtype=float),
            "representative_gene": representative_gene,
            "representative_gene_id": selected_gene_id,
            "effect_role": str(effect.effect_role),
            "effect_role_label": str(effect.effect_role_label),
            "target_state": int(effect.target_state),
            "configured_strength": float(effect.configured_strength),
            "direction": str(effect.direction),
            "selection_rank": int(effect.selection_rank),
            "selection_visualization_score": float(effect.visualization_score),
            "selection_truth_q90_span": float(effect.truth_q90_span),
            "selection_expected_visible_fraction": float(
                effect.expected_visible_fraction
            ),
            "selection_expected_lower_quintile_count": float(
                effect.expected_lower_quintile_count
            ),
            "selection_expected_upper_quintile_count": float(
                effect.expected_upper_quintile_count
            ),
            "selection_poisson_working_contrast_proxy": float(
                effect.poisson_working_contrast_proxy
            ),
            "selection_sqrt_poisson_working_information_proxy": float(
                effect.sqrt_poisson_working_information_proxy
            ),
            "selection_direct_effect_fraction_of_log_count_sd": float(
                effect.direct_effect_fraction_of_log_count_sd
            ),
            "selection_target_high_to_nontarget_q95": float(
                effect.target_high_to_nontarget_q95
            ),
            "is_target_state": (
                simulation.meta["state"].to_numpy(dtype=int) == int(effect.target_state)
            ),
            "direct_log_fold_truth": np.log(latent / base),
            "latent_mean": latent,
            "observed_count": simulation.omics[representative_gene].to_numpy(dtype=float),
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
    seed_summary = effect.groupby(
        ["seed", "basis", "effect_role", "configured_strength"],
        as_index=False,
    )[["latent_strength", "observed_strength"]].median()
    basis_medians = seed_summary.groupby(
        ["basis", "effect_role", "configured_strength"]
    )["observed_strength"].median()
    monotonic_differences = np.concatenate(
        [
            np.diff(
                basis_medians.xs((basis_name, effect_role), level=("basis", "effect_role"))
                .sort_index()
                .to_numpy(dtype=float)
            )
            for basis_name in SPATIAL_BASIS_ORDER
            for effect_role in EFFECT_ROLES
        ]
    )
    configured_values = basis_medians.index.get_level_values("configured_strength").to_numpy(
        dtype=float
    )
    nonzero = effect[effect["is_direct"]].copy()
    sign_recovery = np.mean(
        np.sign(nonzero["configured_strength"].to_numpy())
        == np.sign(nonzero["observed_strength"].to_numpy())
    )
    null_medians = basis_medians.xs(0.0, level="configured_strength")
    null_fpr_by_group = effect.loc[effect["is_null"]].groupby(
        ["basis", "effect_role"]
    )["detected"].mean()
    null_fpr_by_basis = effect.loc[effect["is_null"]].groupby("basis")["detected"].mean()
    powered = effect.loc[
        effect["is_direct"] & (effect["configured_strength"].abs() >= 0.5)
    ].assign(effect_sign=lambda frame: np.sign(frame["configured_strength"]))
    power_by_group = powered.groupby(
        ["basis", "effect_role", "effect_sign"]
    )["detected"].mean()
    maximum_latent_error = float(effect["maximum_latent_spatial_vector_error"].max())
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
            "check": "balanced_effect_and_null_roles",
            "value": float(
                effect.loc[effect["seed"] == min(SEEDS)].groupby(
                    ["basis", "effect_role", "is_null"]
                ).size().min()
            ),
            "criterion": "each basis has 24 nonzero and 6 null genes per role",
            "passed": bool(
                (
                    effect.loc[(effect["seed"] == min(SEEDS)) & effect["is_direct"]]
                    .groupby(["basis", "effect_role"])
                    .size()
                    == 24
                ).all()
                and (
                    effect.loc[(effect["seed"] == min(SEEDS)) & effect["is_null"]]
                    .groupby(["basis", "effect_role"])
                    .size()
                    == 6
                ).all()
            ),
        },
        {
            "check": "fixed_effect_fits_converged",
            "value": float(
                np.mean(effect["full_fit_converged"] & effect["reduced_fit_converged"])
            ),
            "criterion": "equal to 1.0",
            "passed": bool(
                (effect["full_fit_converged"] & effect["reduced_fit_converged"]).all()
            ),
        },
        {
            "check": "fixed_effect_design_full_rank",
            "value": float(
                (effect["fixed_design_rank"] / effect["fixed_design_columns"]).min()
            ),
            "criterion": "rank equals the number of design columns",
            "passed": bool(
                (effect["fixed_design_rank"] == effect["fixed_design_columns"]).all()
            ),
        },
        {
            "check": "full_model_dominates_reduced_model",
            "value": float(effect["full_minus_reduced_log_likelihood"].min()),
            "criterion": "minimum full-minus-reduced log likelihood >= -1e-6",
            "passed": bool(
                effect["full_minus_reduced_log_likelihood"].min() >= -1e-6
            ),
        },
        {
            "check": "lrt_p_values_valid",
            "value": float(effect["p_value"].notna().mean()),
            "criterion": "all p-values are finite and in [0, 1]",
            "passed": bool(
                np.isfinite(effect["p_value"]).all()
                and effect["p_value"].between(0.0, 1.0).all()
            ),
        },
        {
            "check": "latent_spatial_vectors_exact",
            "value": maximum_latent_error,
            "criterion": "maximum absolute error < 1e-10",
            "passed": bool(maximum_latent_error < 1e-10),
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
            "criterion": "at least 0.95 across nonzero gene-level fits",
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
            "check": "null_lrt_type_i_error_upper_bound",
            "value": float(null_fpr_by_group.max()),
            "criterion": (
                "maximum basis-by-role null rejection rate <= 0.125 at "
                "unadjusted per-gene alpha 0.05"
            ),
            "passed": bool(null_fpr_by_group.max() <= 0.125),
        },
        {
            "check": "pooled_null_lrt_type_i_error",
            "value": float(np.max(np.abs(null_fpr_by_basis - DETECTION_ALPHA))),
            "criterion": (
                "each basis-pooled null rejection rate lies between 0.02 and "
                "0.08 at unadjusted per-gene alpha 0.05"
            ),
            "passed": bool(null_fpr_by_basis.between(0.02, 0.08).all()),
        },
        {
            "check": "nonzero_lrt_power_ge_0_5",
            "value": float(power_by_group.min()),
            "criterion": (
                "minimum basis-by-role-by-sign power >= 0.80 for "
                "|coefficient| >= 0.5"
            ),
            "passed": bool(power_by_group.min() >= 0.80),
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
    runtime_version = getattr(ss, "__version__", "unknown")
    if runtime_version != SIMSPACE_EXPECTED_VERSION:
        raise RuntimeError(
            f"Supplementary Figure 3 requires SimSpace {SIMSPACE_EXPECTED_VERSION}; "
            f"imported {runtime_version} from {Path(ss.__file__).resolve()}. "
            "Use PYTHONPATH=../SimSpace to run against the adjacent source checkout."
        )

    effect_rows: list[dict[str, object]] = []
    calibration_rows: list[dict[str, object]] = []
    observation_rows: list[dict[str, object]] = []
    compatibility_rows: list[dict[str, object]] = []
    truth_rows: list[dict[str, object]] = []
    representative: pd.DataFrame | None = None
    representative_effect: pd.Series | None = None
    representative_screen: pd.DataFrame | None = None

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
            effect_rows.extend(fit_spatial_effects(simulation, basis_name).to_dict("records"))

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
        if seed == REPRESENTATIVE_SEED:
            representative_effect, representative_screen = select_representative_effect(clean)
        if representative_effect is None:
            raise RuntimeError("The representative effect must be selected from seed 0 first.")
        truth_rows.extend(gene_truth_rows(clean, seed))
        condition_simulations: dict[str, ss.SimSpace] = {"Clean": clean}

        for efficiency in CAPTURE_LEVELS:
            if efficiency == 1.0:
                simulation = clean
            else:
                simulation = apply_observation_variant(
                    clean,
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
                simulation = apply_observation_variant(
                    clean,
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
                simulation = apply_observation_variant(
                    clean,
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
        combined = apply_observation_variant(
            clean,
            technical_noise=combined_technical,
            dropout=combined_dropout,
        )
        condition_simulations["Combined"] = combined

        clean_median_library = float(
            np.median(clean.omics.to_numpy(dtype=float).sum(axis=1))
        )
        for condition in CONDITION_ORDER:
            metrics = expression_metrics(
                condition_simulations[condition],
                condition,
                seed,
                representative_effect,
            )
            metrics["median_library_ratio"] = (
                float(metrics["median_library_size"]) / clean_median_library
            )
            observation_rows.append(metrics)

        if seed == REPRESENTATIVE_SEED:
            representative = representative_map(combined, representative_effect)
        print(
            f"Completed simulation seed {seed} ({SEEDS.index(seed) + 1}/{len(SEEDS)}).",
            flush=True,
        )

    if representative is None or representative_screen is None:
        raise RuntimeError("The representative map and candidate screen were not generated.")

    frames = {
        "effect": pd.DataFrame(effect_rows),
        "calibration": pd.DataFrame(calibration_rows),
        "observation": pd.DataFrame(observation_rows),
        "compatibility": pd.DataFrame(compatibility_rows),
        "truth": pd.DataFrame(truth_rows),
        "representative": representative,
        "representative_screen": representative_screen,
    }
    checks = acceptance_checks(
        frames["effect"], frames["calibration"], frames["compatibility"]
    )
    frames["checks"] = checks
    if not checks["passed"].all():
        failures = checks.loc[~checks["passed"], ["check", "value", "criterion"]]
        raise RuntimeError(f"Acceptance checks failed:\n{failures.to_string(index=False)}")

    effect_seed_summary = frames["effect"].groupby(
        ["seed", "basis", "effect_scope", "effect_role", "configured_strength"],
        as_index=False,
    )[["latent_strength", "observed_strength"]].median()
    detection_frame = frames["effect"].assign(
        absolute_strength=frames["effect"]["configured_strength"].abs()
    )
    detection_seed_summary = detection_frame.groupby(
        ["seed", "basis", "effect_scope", "effect_role", "absolute_strength"],
        as_index=False,
    )["detected"].mean().rename(columns={"detected": "detection_rate"})
    summary_rows: list[dict[str, object]] = []
    summary_rows.extend(
        summarize_values(
            effect_seed_summary,
            ["basis", "effect_scope", "effect_role", "configured_strength"],
            ["latent_strength", "observed_strength"],
            "direct_spatial_effect",
        )
    )
    summary_rows.extend(
        summarize_values(
            detection_seed_summary,
            ["basis", "effect_scope", "effect_role", "absolute_strength"],
            ["detection_rate"],
            "block_lrt_detection",
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
    write_configuration(frames)
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
    frames["representative_screen"].to_csv(
        DATA_DIR / "representative_candidate_screen.tsv",
        sep="\t",
        index=False,
    )
    frames["checks"].to_csv(DATA_DIR / "acceptance_checks.tsv", sep="\t", index=False)
    frames["summary"].to_csv(DATA_DIR / "summary_metrics.tsv", sep="\t", index=False)


def write_configuration(frames: dict[str, pd.DataFrame]) -> None:
    representative = frames["representative"].iloc[0]
    simspace_import_path = Path(ss.__file__).resolve()
    adjacent_checkout = (REPO_ROOT.parent / "SimSpace").resolve()
    loaded_from_adjacent_checkout = (
        simspace_import_path == adjacent_checkout
        or adjacent_checkout in simspace_import_path.parents
    )
    configuration = {
        "analysis": "Supplementary Figure 3 finite-basis fixed spatial effects and observation layers",
        "analysis_script_sha256": sha256_file(Path(__file__).resolve()),
        "reviewer_comments": ["R2-M1", "R2-m5"],
        "related_comment": "R2-M3",
        "simspace_expected_version": SIMSPACE_EXPECTED_VERSION,
        "simspace_runtime_version": getattr(ss, "__version__", "unknown"),
        "simspace_effective_distribution_version": package_version("simspace"),
        "simspace_distributions_visible_at_run": visible_package_distributions(
            "simspace"
        ),
        "simspace_import_path": str(simspace_import_path),
        "simspace_loaded_from_adjacent_checkout": loaded_from_adjacent_checkout,
        "adjacent_simspace_commit_at_run": adjacent_simspace_commit(),
        "adjacent_simspace_source_dirty_at_run": adjacent_simspace_source_dirty(),
        "adjacent_simspace_source_sha256": {
            str(path.relative_to(adjacent_checkout)): sha256_file(path)
            for path in (
                adjacent_checkout / "simspace" / "omics.py",
                adjacent_checkout / "simspace" / "core.py",
                adjacent_checkout / "simspace" / "__init__.py",
                adjacent_checkout / "tests" / "test_omics.py",
            )
            if path.exists()
        },
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
        "representative_gene_id": int(representative["representative_gene_id"]),
        "representative_gene": str(representative["representative_gene"]),
        "representative_gene_selection": REPRESENTATIVE_GENE_SELECTION,
        "representative_gene_selection_metrics": {
            "rank": int(representative["selection_rank"]),
            "visualization_score": float(
                representative["selection_visualization_score"]
            ),
            "truth_q90_span": float(representative["selection_truth_q90_span"]),
            "expected_visible_fraction": float(
                representative["selection_expected_visible_fraction"]
            ),
            "expected_lower_quintile_count": float(
                representative["selection_expected_lower_quintile_count"]
            ),
            "expected_upper_quintile_count": float(
                representative["selection_expected_upper_quintile_count"]
            ),
            "poisson_working_contrast_proxy": float(
                representative["selection_poisson_working_contrast_proxy"]
            ),
            "sqrt_poisson_working_information_proxy": float(
                representative[
                    "selection_sqrt_poisson_working_information_proxy"
                ]
            ),
            "direct_effect_fraction_of_log_count_sd": float(
                representative["selection_direct_effect_fraction_of_log_count_sd"]
            ),
            "target_high_to_nontarget_q95": float(
                representative["selection_target_high_to_nontarget_q95"]
            ),
        },
        "representative_gene_selection_weights": REPRESENTATIVE_SCORE_WEIGHTS,
        "representative_gene_effect": {
            "basis": "linear",
            "effect_role": str(representative["effect_role"]),
            "effect_role_label": str(representative["effect_role_label"]),
            "target_state": int(representative["target_state"]),
            "configured_strength": float(representative["configured_strength"]),
            "direction": str(representative["direction"]),
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
            "linear_validation_direction": SPATIAL_DIRECTION.tolist(),
            "radial_center": [0.5, 0.5],
            "hotspot_center": [0.30, 0.70],
            "hotspot_bandwidth": 0.18,
            "structure_points": direct_spatial_config("structure_distance")[
                "structure_points"
            ],
            "coordinate_scaling": "unit_range",
            "center_basis": True,
        },
        "spatial_strengths": list(SPATIAL_STRENGTHS),
        "n_effect_genes_per_basis": N_EFFECT_GENES,
        "n_null_validation_genes_per_basis": N_NULL_VALIDATION_GENES,
        "reference_state": REFERENCE_STATE,
        "effect_roles": [
            {
                "role": role,
                "label": EFFECT_ROLE_LABELS[role],
                "nonzero_genes_per_basis": 24,
                "null_genes_per_basis": N_NULLS_PER_ROLE,
            }
            for role in EFFECT_ROLES
        ],
        "nonzero_genes_per_strength_and_role": N_EFFECT_REPLICATES,
        "fixed_effect_model": {
            "link": "Poisson log link",
            "formula": (
                "intercept + non-reference phenotype main effects + B(s) + "
                "non-reference phenotype-by-B(s) interactions"
            ),
            "library_size_factor": 1.0,
            "estimation": "maximum likelihood by Fisher scoring; no oracle baseline offset",
            "test": "target-block likelihood-ratio test against the nested model",
            "test_degrees_of_freedom": "basis dimension (2 for linear; 1 otherwise)",
            "detection_alpha": DETECTION_ALPHA,
        },
        "capture_levels": list(CAPTURE_LEVELS),
        "ambient_rates_counts_per_cell": list(AMBIENT_LEVELS),
        "dropout_levels": [
            {"label": label, "configuration": dropout_configuration}
            for label, dropout_configuration in DROPOUT_LEVELS
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
        "representative_screen": DATA_DIR / "representative_candidate_screen.tsv",
        "checks": DATA_DIR / "acceptance_checks.tsv",
        "summary": DATA_DIR / "summary_metrics.tsv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing archived Supplementary Figure 3 data: " + ", ".join(missing))
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
    representative_gene = str(representative["representative_gene"].iloc[0])
    representative_display = representative_gene.replace("_", " ")
    target_state = int(representative["target_state"].iloc[0])
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
    axes[1].set_title(f"Type {target_state}-specific log-fold truth", fontsize=9)
    colorbar = fig.colorbar(truth_plot, cax=support_axes[1], orientation="horizontal")
    colorbar.set_label("Interaction B(s)ᵀη", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6, length=2)

    target_mask = state == target_state
    axes[2].scatter(
        representative.loc[~target_mask, "col"],
        representative.loc[~target_mask, "row"],
        c="#D3D3D3",
        s=4,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    observed_plot = axes[2].scatter(
        representative.loc[target_mask, "col"],
        representative.loc[target_mask, "row"],
        c=np.log1p(representative.loc[target_mask, "observed_count"]),
        cmap="viridis",
        vmin=0.0,
        s=4,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )
    axes[2].set_title(
        f"Observed {representative_display} in Type {target_state} cells",
        fontsize=9,
    )
    colorbar = fig.colorbar(observed_plot, cax=support_axes[2], orientation="horizontal")
    colorbar.set_label(
        "log(1 + count); target type only, other types masked gray",
        fontsize=6.5,
    )
    colorbar.ax.tick_params(labelsize=6, length=2)

    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    add_panel_label(fig, spec, "A")


def draw_recovery_axis(
    axis: plt.Axes,
    effect: pd.DataFrame,
    effect_scope: str,
    title: str,
) -> None:
    seed_summary = effect.loc[effect["effect_scope"] == effect_scope].groupby(
        ["seed", "basis", "configured_strength"],
        as_index=False,
    )[["latent_strength", "observed_strength"]].median()
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
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_xlabel("Configured coefficient", fontsize=7.5)
    axis.set_title(title, fontsize=8.5)
    axis.tick_params(labelsize=6.5)
    axis.grid(color="#DDDDDD", linewidth=0.45, alpha=0.7)


def draw_panel_b(fig: plt.Figure, spec, effect: pd.DataFrame) -> None:
    grid = spec.subgridspec(
        2,
        2,
        height_ratios=(1.12, 0.88),
        hspace=0.62,
        wspace=0.34,
    )
    overall_axis = fig.add_subplot(grid[0, 0])
    interaction_axis = fig.add_subplot(grid[0, 1])
    detection_axis = fig.add_subplot(grid[1, :])
    draw_recovery_axis(
        overall_axis,
        effect,
        "overall",
        "Base/reference spatial term f₀",
    )
    draw_recovery_axis(
        interaction_axis,
        effect,
        "cell_type_specific",
        "Cell-type-specific deviations fₖ",
    )
    overall_axis.set_ylabel("Recovered coefficient", fontsize=7.5)
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
    interaction_axis.legend(
        handles=basis_handles,
        frameon=False,
        fontsize=5.6,
        loc="upper left",
        ncol=2,
        columnspacing=0.5,
        handletextpad=0.25,
    )

    detection = effect.assign(
        absolute_strength=effect["configured_strength"].abs()
    ).groupby(
        ["seed", "effect_scope", "absolute_strength"],
        as_index=False,
    )["detected"].mean().rename(columns={"detected": "detection_rate"})
    detection_styles = {
        "overall": ("#222222", "o", "Base/reference f₀"),
        "cell_type_specific": ("#7A5195", "^", "Cell-type-specific fₖ"),
    }
    rng = np.random.default_rng(418)
    for scope, (color, marker, label) in detection_styles.items():
        scope_data = detection[detection["effect_scope"] == scope]
        x_values = np.sort(scope_data["absolute_strength"].unique().astype(float))
        medians: list[float] = []
        for x_value in x_values:
            values = scope_data.loc[
                scope_data["absolute_strength"] == x_value,
                "detection_rate",
            ].to_numpy(dtype=float)
            jitter = rng.normal(0.0, 0.008, size=len(values))
            detection_axis.scatter(
                np.full(len(values), x_value) + jitter,
                values,
                s=8,
                marker=marker,
                color=color,
                alpha=0.22,
                edgecolors="none",
                rasterized=True,
            )
            median = float(np.median(values))
            medians.append(median)
            detection_axis.errorbar(
                x_value,
                median,
                yerr=[[median - values.min()], [values.max() - median]],
                fmt=marker,
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.5,
                markersize=4.5,
                linewidth=0.9,
                capsize=1.5,
                zorder=5,
            )
        detection_axis.plot(x_values, medians, color=color, linewidth=1.0, label=label)
    detection_axis.axhline(
        DETECTION_ALPHA,
        color="#777777",
        linestyle="--",
        linewidth=0.8,
        label="Nominal α",
    )
    detection_axis.set_xlim(-0.055, 1.055)
    detection_axis.set_ylim(-0.025, 1.035)
    detection_axis.set_xticks((0.0, 0.25, 0.5, 0.75, 1.0))
    detection_axis.set_xlabel(
        "|Configured coefficient| (0 = null rejection rate)",
        fontsize=7.5,
    )
    detection_axis.set_ylabel("Detection rate per seed", fontsize=7.5)
    detection_axis.set_title(
        "Target-block likelihood-ratio test (df = basis dimension)",
        fontsize=8.5,
        loc="left",
    )
    detection_axis.tick_params(labelsize=6.5)
    detection_axis.grid(color="#E3E3E3", linewidth=0.45)
    detection_axis.legend(
        frameon=False,
        fontsize=6.0,
        loc="lower right",
        ncol=3,
        columnspacing=0.7,
        handletextpad=0.3,
    )
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
    png_path = SCRIPT_DIR / "SFig3.png"
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    shutil.copy2(png_path, EXAMPLE_DIR / png_path.name)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--render-only",
        action="store_true",
        help="Render from archived TSV files without rerunning the simulations.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    frames = load_data() if arguments.render_only else generate_data()
    render_figure(frames)
    print(f"Wrote {SCRIPT_DIR / 'SFig3.png'}")
    if not arguments.render_only:
        print(f"All {len(frames['checks'])} acceptance checks passed.")


if __name__ == "__main__":
    main()
