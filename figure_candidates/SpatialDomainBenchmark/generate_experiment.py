#!/usr/bin/env python3
"""Generate the frozen 2 x 3 x 3 spatial-domain benchmark design."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from generate_pilot import (
    N_CELL_TYPES,
    N_DOMAINS,
    N_GENES,
    N_LOCATIONS,
    generate_dataset,
)
from pilot_common import sha256_file


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_ROOT = SCRIPT_DIR / "experiment_data"
TOPOLOGIES = ("curved_layers", "irregular_mrf")
DIFFICULTY_FOLD_CHANGES = {
    "hard": 2.0,
    "moderate": 3.0,
    "easy": 4.0,
}
LAYOUT_SEEDS = (0, 1, 2)
METHODS = ("graphst", "stagate", "banksy", "spagcn", "spclue", "stlearn")
NEGATIVE_CONTROLS = ("expression_pca_kmeans",)
PAIRED_FILES = ("coordinates.tsv", "truth.tsv", "genes.tsv", "gene_truth.tsv")


def dataset_name(topology: str, seed: int, difficulty: str) -> str:
    return f"{topology}_seed{seed}_{difficulty}"


def expected_dataset_names(
    topologies: tuple[str, ...],
    difficulties: tuple[str, ...],
    layout_seeds: tuple[int, ...],
) -> tuple[str, ...]:
    return tuple(
        dataset_name(topology, seed, difficulty)
        for topology in topologies
        for seed in layout_seeds
        for difficulty in difficulties
    )


def validate_dataset(
    dataset_dir: Path,
    topology: str,
    difficulty: str,
    seed: int,
) -> dict[str, object]:
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing dataset manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    expected = {
        "dataset": dataset_dir.name,
        "topology": topology,
        "difficulty": difficulty,
        "layout_seed": seed,
        "n_locations": N_LOCATIONS,
        "n_genes": N_GENES,
        "n_domains": N_DOMAINS,
        "n_cell_types": N_CELL_TYPES,
        "domain_fold_change": DIFFICULTY_FOLD_CHANGES[difficulty],
    }
    for field, value in expected.items():
        if manifest.get(field) != value:
            raise ValueError(
                f"{manifest_path}: expected {field}={value!r}, "
                f"found {manifest.get(field)!r}"
            )
    checksums = manifest.get("checksums", {})
    for filename, expected_checksum in checksums.items():
        path = dataset_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Missing generated file: {path}")
        observed_checksum = sha256_file(path)
        if observed_checksum != expected_checksum:
            raise ValueError(f"Checksum mismatch for {path}")
    return manifest


def validate_paired_difficulties(
    output_root: Path,
    topologies: tuple[str, ...],
    difficulties: tuple[str, ...],
    layout_seeds: tuple[int, ...],
) -> None:
    """Ensure only the domain-expression signal changes within a paired trio."""
    for topology in topologies:
        for seed in layout_seeds:
            directories = [
                output_root / dataset_name(topology, seed, difficulty)
                for difficulty in difficulties
            ]
            for filename in PAIRED_FILES:
                checksums = {sha256_file(directory / filename) for directory in directories}
                if len(checksums) != 1:
                    raise ValueError(
                        f"Paired file {filename} differs across difficulty levels "
                        f"for {topology}, layout seed {seed}"
                    )
            for filename in ("latent_mean.npy.gz", "counts.mtx.gz"):
                checksums = {
                    sha256_file(directory / filename) for directory in directories
                }
                if len(difficulties) > 1 and len(checksums) != len(difficulties):
                    raise ValueError(
                        f"{filename} is not distinct across all difficulty levels "
                        f"for {topology}, layout seed {seed}"
                    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate the reference-free production benchmark datasets."
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--topologies", nargs="+", choices=TOPOLOGIES, default=TOPOLOGIES
    )
    parser.add_argument(
        "--difficulties",
        nargs="+",
        choices=tuple(DIFFICULTY_FOLD_CHANGES),
        default=tuple(DIFFICULTY_FOLD_CHANGES),
    )
    parser.add_argument(
        "--layout-seeds", nargs="+", type=int, choices=LAYOUT_SEEDS, default=LAYOUT_SEEDS
    )
    args = parser.parse_args()

    topologies = tuple(dict.fromkeys(args.topologies))
    difficulties = tuple(dict.fromkeys(args.difficulties))
    layout_seeds = tuple(dict.fromkeys(args.layout_seeds))
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    dataset_records: list[dict[str, object]] = []
    for topology in topologies:
        for seed in layout_seeds:
            for difficulty in difficulties:
                dataset_dir = generate_dataset(
                    output_root,
                    topology=topology,
                    seed=seed,
                    difficulty=difficulty,
                    domain_fold_change=DIFFICULTY_FOLD_CHANGES[difficulty],
                )
                manifest = validate_dataset(
                    dataset_dir,
                    topology=topology,
                    difficulty=difficulty,
                    seed=seed,
                )
                dataset_records.append(
                    {
                        "dataset": dataset_dir.name,
                        "topology": topology,
                        "difficulty": difficulty,
                        "layout_seed": seed,
                        "domain_fold_change": DIFFICULTY_FOLD_CHANGES[difficulty],
                        "manifest_sha256": sha256_file(dataset_dir / "manifest.json"),
                        "counts_sha256": manifest["checksums"]["counts.mtx.gz"],
                    }
                )
                print(f"Wrote {dataset_dir}", flush=True)

    validate_paired_difficulties(
        output_root, topologies, difficulties, layout_seeds
    )
    expected_names = expected_dataset_names(topologies, difficulties, layout_seeds)
    observed_names = tuple(record["dataset"] for record in dataset_records)
    if observed_names != expected_names:
        raise RuntimeError("Generated dataset order does not match the frozen design")

    root_manifest = {
        "schema_version": "1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "design": "2 pattern types x 3 difficulty levels x 3 layout seeds",
        "topologies": list(topologies),
        "difficulties": {
            difficulty: DIFFICULTY_FOLD_CHANGES[difficulty]
            for difficulty in difficulties
        },
        "layout_seeds": list(layout_seeds),
        "n_datasets": len(dataset_records),
        "named_methods": list(METHODS),
        "n_named_method_runs": len(dataset_records) * len(METHODS),
        "negative_controls": list(NEGATIVE_CONTROLS),
        "n_negative_control_runs": len(dataset_records) * len(NEGATIVE_CONTROLS),
        "n_total_evaluated_runs": len(dataset_records)
        * (len(METHODS) + len(NEGATIVE_CONTROLS)),
        "paired_across_difficulty": {
            "within": "topology and layout_seed",
            "byte_identical_files": list(PAIRED_FILES),
            "changing_factor": "nominal domain fold change",
        },
        "datasets": dataset_records,
    }
    (output_root / "experiment_manifest.json").write_text(
        json.dumps(root_manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Validated {len(dataset_records)} datasets and "
        f"{len(dataset_records) * len(METHODS)} planned named-method runs."
    )


if __name__ == "__main__":
    main()
