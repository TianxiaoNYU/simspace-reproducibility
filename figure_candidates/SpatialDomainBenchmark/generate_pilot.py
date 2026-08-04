#!/usr/bin/env python3
"""Generate truth-separated reference-free pilot datasets."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE = SCRIPT_DIR.parents[3]
LOCAL_SIMSPACE = WORKSPACE / "SimSpace"
if str(LOCAL_SIMSPACE) not in sys.path:
    sys.path.insert(0, str(LOCAL_SIMSPACE))

import simspace as ss  # noqa: E402

from pilot_common import sha256_file  # noqa: E402


N_LOCATIONS = 3000
N_GENES = 2000
N_DOMAINS = 4
N_CELL_TYPES = 10
DEFAULT_PILOT_DIFFICULTY = "intermediate"
DEFAULT_DOMAIN_FOLD_CHANGE = 3.0
ENRICHED_TYPES_PER_DOMAIN = (3, 2, 3, 2)
ENRICHED_MASS = 0.40
SHARED_MASS = 0.60
PHENOTYPE_FOLD_CHANGE = 4.0
DOMAIN_PROGRAM_GENES_PER_DOMAIN = 50
DOMAIN_PROGRAM_START = N_GENES - N_DOMAINS * DOMAIN_PROGRAM_GENES_PER_DOMAIN
DOMAIN_STRENGTH_TEMPLATE = (0.80, 0.95, 1.05, 1.20)
DOMAIN_GENE_WEIGHT_RANGE = (0.60, 1.40)
DOMAIN_ACTIVITY_RANGE = (0.50, 1.00)
DISPERSION = 8.0
LIBRARY_SIGMA = 0.25
IRREGULAR_GRID_SHAPE = (68, 68)
IRREGULAR_COORDINATE_JITTER = 0.55
NICHE_MRF_SEED_BASE = 7008
NICHE_MRF_PHI = 5.25
NICHE_MRF_SWEEPS = 40
NICHE_MRF_NEIGHBORHOOD_RADIUS = 2
NICHE_MRF_DIAGONAL = 1.0
NICHE_MRF_OFF_DIAGONAL = -0.12


def smooth_noise(shape: tuple[int, int], seed: int, sigma: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    values = gaussian_filter(rng.normal(size=shape), sigma=sigma, mode="reflect")
    return values / values.std()


def curved_layers(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Create a reference-free curved tissue mask and four radial domains."""
    rows, cols = 104, 118
    row, col = np.indices((rows, cols))
    x = col / (cols - 1)
    y = row / (rows - 1)
    dx = x - 1.15
    dy = (y - 1.18) / 1.05
    radius = np.sqrt(dx**2 + dy**2)
    angle = np.degrees(np.arctan2(dy, dx))
    angular_position = np.clip((-82.0 - angle) / 87.0, 0.0, 1.0)
    outer = 1.14 + 0.18 * angular_position
    outer += 0.010 * smooth_noise((rows, cols), 419 + seed, 5.0)
    tissue = (
        (angle >= -169.0)
        & (angle <= -82.0)
        & (radius >= 0.44)
        & (radius <= outer)
    )
    depth = np.clip((outer - radius) / (outer - 0.44), 0.0, 1.0)
    interface = 0.018 * smooth_noise((rows, cols), 1420 + seed, 7.0)
    cuts = np.sort(
        np.stack(
            [
                np.full_like(depth, 0.29) + interface,
                np.full_like(depth, 0.51) + 0.65 * interface,
                np.full_like(depth, 0.69) + 0.45 * interface,
            ]
        ),
        axis=0,
    )
    domains = np.zeros((rows, cols), dtype=int)
    domains[depth >= cuts[0]] = 1
    domains[depth >= cuts[1]] = 2
    domains[depth >= cuts[2]] = 3
    domains[~tissue] = 0
    return tissue, domains


def irregular_mrf(seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Create four irregular domains using SimSpace's niche-level MRF."""
    shape = IRREGULAR_GRID_SHAPE
    placeholder_theta = [np.eye(N_CELL_TYPES)] * N_DOMAINS
    simulation = ss.SimSpace(
        shape=shape,
        num_states=N_CELL_TYPES,
        theta=placeholder_theta,
        num_iterations=1,
        phi=NICHE_MRF_PHI,
        random_seed=NICHE_MRF_SEED_BASE + seed,
    )
    simulation.initialize()
    niche_theta = np.full(
        (N_DOMAINS, N_DOMAINS), NICHE_MRF_OFF_DIAGONAL, dtype=float
    )
    np.fill_diagonal(niche_theta, NICHE_MRF_DIAGONAL)
    simulation.create_niche(
        num_niches=N_DOMAINS,
        n_iter=NICHE_MRF_SWEEPS,
        theta_niche=niche_theta,
        neighborhood=ss.spatial.generate_offsets(
            NICHE_MRF_NEIGHBORHOOD_RADIUS, method="manhattan"
        ),
    )
    return np.ones(shape, dtype=bool), simulation.niche.copy()


def home_type_sets(seed: int) -> list[np.ndarray]:
    permutation = np.roll(np.arange(N_CELL_TYPES), seed % N_CELL_TYPES)
    groups: list[np.ndarray] = []
    start = 0
    for size in ENRICHED_TYPES_PER_DOMAIN:
        groups.append(permutation[start : start + size])
        start += size
    return groups


def cell_type_probabilities(seed: int) -> np.ndarray:
    probabilities = np.full(
        (N_DOMAINS, N_CELL_TYPES), SHARED_MASS / N_CELL_TYPES, dtype=float
    )
    for domain, home_types in enumerate(home_type_sets(seed)):
        probabilities[domain, home_types] += ENRICHED_MASS / len(home_types)
    return probabilities


def domain_strengths(seed: int) -> np.ndarray:
    """Rotate unequal domain strengths so no label is systematically easiest."""
    return np.roll(np.asarray(DOMAIN_STRENGTH_TEMPLATE, dtype=float), seed % N_DOMAINS)


def domain_gene_effect_weights(seed: int) -> np.ndarray:
    """Generate fixed heterogeneous log-effect weights for domain genes."""
    rng = np.random.default_rng(11000 + seed)
    low, high = DOMAIN_GENE_WEIGHT_RANGE
    weights = rng.beta(
        2.0,
        2.0,
        size=(N_DOMAINS, DOMAIN_PROGRAM_GENES_PER_DOMAIN),
    )
    return low + (high - low) * weights


def domain_location_activity(seed: int) -> np.ndarray:
    """Generate fixed partial activation for each simulated location."""
    rng = np.random.default_rng(12000 + seed)
    low, high = DOMAIN_ACTIVITY_RANGE
    activity = rng.beta(5.0, 2.0, size=N_LOCATIONS)
    return low + (high - low) * activity


def simulate_locations(
    topology: str, seed: int
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    if topology == "curved_layers":
        tissue, domains = curved_layers(seed)
    elif topology == "irregular_mrf":
        tissue, domains = irregular_mrf(seed)
    else:
        raise ValueError(f"Unknown topology: {topology}")

    available = np.flatnonzero(tissue)
    if len(available) < N_LOCATIONS:
        raise RuntimeError(f"{topology} exposes only {len(available)} locations")
    rng = np.random.default_rng(8000 + seed)
    retained = np.sort(rng.choice(available, size=N_LOCATIONS, replace=False))
    retained_mask = np.zeros(tissue.shape, dtype=bool)
    retained_mask.flat[retained] = True

    domain_counts = np.bincount(domains.flat[retained], minlength=N_DOMAINS)
    if np.any(domain_counts < 0.10 * N_LOCATIONS):
        raise RuntimeError(f"Degenerate domain sizes for {topology}: {domain_counts}")

    probabilities = cell_type_probabilities(seed)
    grid = np.full(tissue.shape, -1, dtype=int)
    for domain in range(N_DOMAINS):
        indices = np.flatnonzero(retained_mask & (domains == domain))
        grid.flat[indices] = rng.choice(
            N_CELL_TYPES, size=len(indices), p=probabilities[domain]
        )

    phenotype_theta = np.full((N_CELL_TYPES, N_CELL_TYPES), 0.02, dtype=float)
    np.fill_diagonal(phenotype_theta, 0.35)
    simulation = ss.SimSpace(
        shape=tissue.shape,
        num_states=N_CELL_TYPES,
        theta=[phenotype_theta.copy() for _ in range(N_DOMAINS)],
        num_iterations=1,
        phi=1.4,
        neighborhood=ss.spatial.generate_offsets(1, method="manhattan"),
        random_seed=9000 + seed,
    )
    simulation.grid = grid
    simulation.niche = domains
    simulation.gibbs_sampler()
    jitter = IRREGULAR_COORDINATE_JITTER if topology == "irregular_mrf" else 0.15
    simulation.perturbation(step=jitter)

    coordinates = simulation.meta[["col", "row"]].to_numpy(dtype=float)
    nearest = cKDTree(coordinates).query(coordinates, k=2)[0][:, 1]
    coordinate_scale = float(np.median(nearest))
    coordinates /= coordinate_scale
    location_ids = np.asarray([f"location_{index:04d}" for index in range(N_LOCATIONS)])
    coordinate_frame = pd.DataFrame(
        {"location_id": location_ids, "x": coordinates[:, 0], "y": coordinates[:, 1]}
    )
    domain_truth = simulation.meta["niche"].astype(int).to_numpy()
    cell_type_truth = simulation.meta["state"].astype(int).to_numpy()
    return coordinate_frame, domain_truth, cell_type_truth


def simulate_counts(
    domain_truth: np.ndarray,
    cell_type_truth: np.ndarray,
    seed: int,
    domain_fold_change: float = DEFAULT_DOMAIN_FOLD_CHANGE,
    return_latent_mean: bool = False,
) -> (
    tuple[sparse.csr_matrix, pd.DataFrame]
    | tuple[sparse.csr_matrix, pd.DataFrame, np.ndarray]
):
    rng = np.random.default_rng(10000 + seed)
    base_mean = rng.gamma(shape=2.0, scale=0.60, size=N_GENES)
    base_mean = np.clip(base_mean, 0.05, None)
    library_factor = rng.lognormal(
        mean=-(LIBRARY_SIGMA**2) / 2,
        sigma=LIBRARY_SIGMA,
        size=N_LOCATIONS,
    )
    mean = library_factor[:, None] * base_mean[None, :]

    gene_class = np.full(N_GENES, "background", dtype=object)
    target = np.full(N_GENES, "none", dtype=object)
    domain_effect_weight = np.full(N_GENES, np.nan, dtype=float)
    target_domain_strength = np.full(N_GENES, np.nan, dtype=float)
    for cell_type in range(N_CELL_TYPES):
        start = 1000 + 60 * cell_type
        stop = start + 60
        gene_class[start:stop] = "phenotype_program"
        target[start:stop] = f"cell_type_{cell_type}"
        rows = cell_type_truth == cell_type
        mean[np.ix_(rows, np.arange(start, stop))] *= PHENOTYPE_FOLD_CHANGE
    gene_weights = domain_gene_effect_weights(seed)
    location_activity = domain_location_activity(seed)
    strengths = domain_strengths(seed)
    for domain in range(N_DOMAINS):
        start = DOMAIN_PROGRAM_START + DOMAIN_PROGRAM_GENES_PER_DOMAIN * domain
        stop = start + DOMAIN_PROGRAM_GENES_PER_DOMAIN
        gene_class[start:stop] = "domain_program"
        target[start:stop] = f"domain_{domain}"
        domain_effect_weight[start:stop] = gene_weights[domain]
        target_domain_strength[start:stop] = strengths[domain]
        rows = domain_truth == domain
        log_multiplier = (
            np.log(domain_fold_change)
            * strengths[domain]
            * location_activity[rows, None]
            * gene_weights[domain][None, :]
        )
        mean[np.ix_(rows, np.arange(start, stop))] *= np.exp(log_multiplier)

    latent_rate = rng.gamma(shape=DISPERSION, scale=mean / DISPERSION)
    counts = rng.poisson(latent_rate).astype(np.int32)
    genes = pd.DataFrame(
        {
            "gene_id": [f"gene_{index:04d}" for index in range(N_GENES)],
            "gene_class": gene_class,
            "target": target,
            "baseline_mean": base_mean,
            "domain_effect_weight": domain_effect_weight,
            "target_domain_strength": target_domain_strength,
        }
    )
    result = sparse.csr_matrix(counts), genes
    if return_latent_mean:
        return (*result, mean.astype(np.float32))
    return result


def write_matrix_market_gzip(matrix: sparse.csr_matrix, path: Path) -> None:
    temporary = path.with_suffix("")
    mmwrite(temporary, matrix, field="integer", symmetry="general")
    with temporary.open("rb") as source, path.open("wb") as compressed:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=compressed, mtime=0
        ) as destination:
            shutil.copyfileobj(source, destination)
    temporary.unlink()


def write_numpy_gzip(array: np.ndarray, path: Path) -> None:
    """Write an auditable dense array with deterministic gzip metadata."""
    with path.open("wb") as compressed:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=compressed, mtime=0
        ) as destination:
            np.save(destination, array, allow_pickle=False)


def generate_dataset(
    root: Path,
    topology: str,
    seed: int,
    difficulty: str = DEFAULT_PILOT_DIFFICULTY,
    domain_fold_change: float = DEFAULT_DOMAIN_FOLD_CHANGE,
) -> Path:
    dataset_name = f"{topology}_seed{seed}_{difficulty}"
    dataset_dir = root / dataset_name
    dataset_dir.mkdir(parents=True, exist_ok=True)
    coordinates, domains, cell_types = simulate_locations(topology, seed)
    counts, gene_truth, latent_mean = simulate_counts(
        domains,
        cell_types,
        seed,
        domain_fold_change=domain_fold_change,
        return_latent_mean=True,
    )

    write_matrix_market_gzip(counts, dataset_dir / "counts.mtx.gz")
    write_numpy_gzip(latent_mean, dataset_dir / "latent_mean.npy.gz")
    gene_truth[["gene_id"]].to_csv(dataset_dir / "genes.tsv", sep="\t", index=False)
    coordinates.to_csv(dataset_dir / "coordinates.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "location_id": coordinates["location_id"],
            "domain_truth": domains,
            "cell_type_truth": cell_types,
            "domain_activity": domain_location_activity(seed),
        }
    ).to_csv(dataset_dir / "truth.tsv", sep="\t", index=False)
    gene_truth.to_csv(dataset_dir / "gene_truth.tsv", sep="\t", index=False)

    files = [
        "counts.mtx.gz",
        "latent_mean.npy.gz",
        "genes.tsv",
        "coordinates.tsv",
        "truth.tsv",
        "gene_truth.tsv",
    ]
    manifest = {
        "dataset": dataset_name,
        "topology": topology,
        "difficulty": difficulty,
        "layout_seed": seed,
        "n_locations": N_LOCATIONS,
        "n_genes": N_GENES,
        "n_domains": N_DOMAINS,
        "n_cell_types": N_CELL_TYPES,
        "latent_mean_format": "gzip-compressed NumPy float32 array, locations_by_genes",
        "domain_fold_change": domain_fold_change,
        "domain_program_genes_per_domain": DOMAIN_PROGRAM_GENES_PER_DOMAIN,
        "domain_program_total_genes": N_DOMAINS * DOMAIN_PROGRAM_GENES_PER_DOMAIN,
        "domain_effect_model": "heterogeneous_log_scale",
        "domain_strengths": domain_strengths(seed).tolist(),
        "domain_gene_weight_range": list(DOMAIN_GENE_WEIGHT_RANGE),
        "domain_gene_weight_distribution": "scaled_beta_2_2",
        "domain_activity_range": list(DOMAIN_ACTIVITY_RANGE),
        "domain_activity_distribution": "scaled_beta_5_2",
        "phenotype_fold_change": PHENOTYPE_FOLD_CHANGE,
        "enriched_types_per_domain": list(ENRICHED_TYPES_PER_DOMAIN),
        "enriched_probability_mass": ENRICHED_MASS,
        "shared_probability_mass": SHARED_MASS,
        "grid_shape": (
            list(IRREGULAR_GRID_SHAPE) if topology == "irregular_mrf" else [104, 118]
        ),
        "retention_fraction": (
            float(N_LOCATIONS / np.prod(IRREGULAR_GRID_SHAPE))
            if topology == "irregular_mrf"
            else None
        ),
        "coordinate_jitter": (
            IRREGULAR_COORDINATE_JITTER if topology == "irregular_mrf" else 0.15
        ),
        "niche_mrf": (
            {
                "seed": NICHE_MRF_SEED_BASE + seed,
                "phi": NICHE_MRF_PHI,
                "sweeps": NICHE_MRF_SWEEPS,
                "neighborhood": "manhattan",
                "neighborhood_radius": NICHE_MRF_NEIGHBORHOOD_RADIUS,
                "theta_diagonal": NICHE_MRF_DIAGONAL,
                "theta_off_diagonal": NICHE_MRF_OFF_DIAGONAL,
            }
            if topology == "irregular_mrf"
            else None
        ),
        "domain_counts": np.bincount(domains, minlength=N_DOMAINS).tolist(),
        "cell_type_counts": np.bincount(cell_types, minlength=N_CELL_TYPES).tolist(),
        "checksums": {name: sha256_file(dataset_dir / name) for name in files},
    }
    (dataset_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return dataset_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=SCRIPT_DIR / "pilot_data")
    parser.add_argument(
        "--topologies",
        nargs="+",
        choices=("curved_layers", "irregular_mrf"),
        default=("irregular_mrf",),
    )
    parser.add_argument(
        "--difficulty-label",
        default=DEFAULT_PILOT_DIFFICULTY,
        help="Label appended to the dataset name and recorded in the manifest.",
    )
    parser.add_argument(
        "--domain-fold-change",
        type=float,
        default=DEFAULT_DOMAIN_FOLD_CHANGE,
    )
    args = parser.parse_args()
    if not args.difficulty_label.replace("_", "").isalnum():
        raise ValueError("difficulty label must contain only letters, digits, or underscores")
    if args.domain_fold_change <= 1.0:
        raise ValueError("domain fold change must be greater than 1")
    args.output_root.mkdir(parents=True, exist_ok=True)
    for topology in args.topologies:
        output = generate_dataset(
            args.output_root,
            topology,
            seed=0,
            difficulty=args.difficulty_label,
            domain_fold_change=args.domain_fold_change,
        )
        print(f"Wrote {output}")


if __name__ == "__main__":
    main()
