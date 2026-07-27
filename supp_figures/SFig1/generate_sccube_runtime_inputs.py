#!/usr/bin/env python3
"""Generate deterministic 2,200-feature scCube runtime inputs.

The generated matrices are synthetic dimensionality stress-test inputs for the
spatial-pattern runtime benchmark. They are not intended to represent a
biologically realistic 2,200-gene expression panel.

For each source gene, the original profile is retained exactly and nine
nonnegative variants are created by adding deterministic Gaussian jitter with
standard deviation max(0.05 * source_gene_sd, 0.001). Outputs are written
atomically and accompanied by a checksum manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import os
import re
import shutil
from pathlib import Path
from typing import Iterator

import numpy as np


DEFAULT_SOURCE_GENES = 220
DEFAULT_COPIES_PER_GENE = 10
DEFAULT_BASE_SEED = 2200
DEFAULT_RELATIVE_NOISE_SD = 0.05
DEFAULT_MINIMUM_NOISE_SD = 1e-3
GRID_PATTERN = re.compile(r"sccube_data_(\d+)\.csv$")


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Generate or verify deterministic expanded scCube inputs for "
            "runtime benchmarking."
        )
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=script_dir / "data",
        help="Directory containing the tracked 220-gene input matrices.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "generated_data",
        help="Directory for generated matrices and copied metadata.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=script_dir / "generated_data_manifest.tsv",
        help="Tracked checksum and parameter manifest.",
    )
    parser.add_argument(
        "--source-genes",
        type=int,
        default=DEFAULT_SOURCE_GENES,
        help="Expected number of source gene rows.",
    )
    parser.add_argument(
        "--copies-per-gene",
        type=int,
        default=DEFAULT_COPIES_PER_GENE,
        help="Total output rows per source gene, including the original.",
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=DEFAULT_BASE_SEED,
        help="Base seed; each file uses base_seed + grid_size.",
    )
    parser.add_argument(
        "--relative-noise-sd",
        type=float,
        default=DEFAULT_RELATIVE_NOISE_SD,
        help="Jitter SD as a fraction of each source gene's across-cell SD.",
    )
    parser.add_argument(
        "--minimum-noise-sd",
        type=float,
        default=DEFAULT_MINIMUM_NOISE_SD,
        help="Minimum jitter SD for low-variance profiles.",
    )
    parser.add_argument(
        "--verify-existing",
        action="store_true",
        help=(
            "Regenerate expected rows in memory and compare them exactly with "
            "the existing output files without rewriting the matrices."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing generated files atomically.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.source_genes <= 0:
        raise ValueError("--source-genes must be positive")
    if args.copies_per_gene < 2:
        raise ValueError("--copies-per-gene must be at least 2")
    if args.relative_noise_sd <= 0:
        raise ValueError("--relative-noise-sd must be positive")
    if args.minimum_noise_sd <= 0:
        raise ValueError("--minimum-noise-sd must be positive")
    if args.verify_existing and args.overwrite:
        raise ValueError("--verify-existing and --overwrite cannot be combined")


def grid_size(path: Path) -> int:
    match = GRID_PATTERN.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unexpected scCube filename: {path.name}")
    return int(match.group(1))


def discover_sources(source_dir: Path) -> list[Path]:
    files = sorted(source_dir.glob("sccube_data_*.csv"), key=grid_size)
    if not files:
        raise FileNotFoundError(
            f"No sccube_data_<grid>.csv files found under {source_dir}"
        )
    return files


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def format_vector(values: np.ndarray) -> str:
    return ",".join(format(value, ".8g") for value in values)


def expected_lines(
    source_path: Path,
    *,
    source_genes: int,
    copies_per_gene: int,
    base_seed: int,
    relative_noise_sd: float,
    minimum_noise_sd: float,
) -> Iterator[str]:
    grid = grid_size(source_path)
    expected_cells = grid * grid
    rng = np.random.default_rng(base_seed + grid)
    seen_names: set[str] = set()
    source_row_count = 0

    with source_path.open("r", encoding="utf-8", newline="") as source:
        header = source.readline()
        if not header:
            raise ValueError(f"{source_path.name}: empty source file")
        header_columns = next(csv.reader([header]))
        if len(header_columns) != expected_cells + 1:
            raise ValueError(
                f"{source_path.name}: expected {expected_cells} cell columns, "
                f"found {len(header_columns) - 1}"
            )
        yield header if header.endswith("\n") else header + "\n"

        for line_number, raw_line in enumerate(source, start=2):
            stripped = raw_line.rstrip("\r\n")
            gene_name, separator, numeric_text = stripped.partition(",")
            if not separator or not gene_name:
                raise ValueError(
                    f"{source_path.name}:{line_number}: malformed source row"
                )
            if gene_name in seen_names:
                raise ValueError(
                    f"{source_path.name}:{line_number}: duplicate gene name "
                    f"{gene_name}"
                )

            values = np.fromstring(numeric_text, dtype=np.float64, sep=",")
            if values.size != expected_cells:
                raise ValueError(
                    f"{source_path.name}:{gene_name}: expected {expected_cells} "
                    f"values, found {values.size}"
                )
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(
                    f"{source_path.name}:{gene_name}: source values must be "
                    "finite and nonnegative"
                )

            yield stripped + "\n"
            seen_names.add(gene_name)
            source_row_count += 1

            gene_sd = float(values.std())
            noise_sd = max(relative_noise_sd * gene_sd, minimum_noise_sd)
            group_hashes = {hashlib.sha256(values.tobytes()).hexdigest()}

            for copy_index in range(1, copies_per_gene):
                variant_name = f"{gene_name}__aug{copy_index:02d}"
                if variant_name in seen_names:
                    raise ValueError(
                        f"{source_path.name}: generated duplicate name "
                        f"{variant_name}"
                    )

                variant = np.maximum(
                    values + rng.normal(0.0, noise_sd, size=values.size),
                    0.0,
                )
                if not np.isfinite(variant).all() or (variant < 0).any():
                    raise ValueError(
                        f"{source_path.name}:{variant_name}: invalid generated "
                        "value"
                    )
                variant_hash = hashlib.sha256(variant.tobytes()).hexdigest()
                if variant_hash in group_hashes:
                    raise ValueError(
                        f"{source_path.name}:{variant_name}: duplicate within "
                        "augmented gene group"
                    )
                group_hashes.add(variant_hash)
                seen_names.add(variant_name)
                yield variant_name + "," + format_vector(variant) + "\n"

    if source_row_count != source_genes:
        raise ValueError(
            f"{source_path.name}: expected {source_genes} source genes, "
            f"found {source_row_count}"
        )
    expected_names = source_genes * copies_per_gene
    if len(seen_names) != expected_names:
        raise ValueError(
            f"{source_path.name}: expected {expected_names} unique output "
            f"names, constructed {len(seen_names)}"
        )


def generate_data_file(
    source_path: Path,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"{output_path} already exists; use --verify-existing to check it "
            "or --overwrite to replace it"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".generating.tmp")

    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as output:
            output.writelines(
                expected_lines(
                    source_path,
                    source_genes=args.source_genes,
                    copies_per_gene=args.copies_per_gene,
                    base_seed=args.base_seed,
                    relative_noise_sd=args.relative_noise_sd,
                    minimum_noise_sd=args.minimum_noise_sd,
                )
            )
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def verify_data_file(
    source_path: Path,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    if not output_path.exists():
        raise FileNotFoundError(f"Missing generated file: {output_path}")
    missing = object()
    expected = expected_lines(
        source_path,
        source_genes=args.source_genes,
        copies_per_gene=args.copies_per_gene,
        base_seed=args.base_seed,
        relative_noise_sd=args.relative_noise_sd,
        minimum_noise_sd=args.minimum_noise_sd,
    )
    with output_path.open("r", encoding="utf-8", newline="") as observed:
        for line_number, (expected_line, observed_line) in enumerate(
            itertools.zip_longest(expected, observed, fillvalue=missing),
            start=1,
        ):
            if expected_line is missing:
                raise ValueError(
                    f"{output_path.name}:{line_number}: unexpected extra row"
                )
            if observed_line is missing:
                raise ValueError(
                    f"{output_path.name}:{line_number}: generated file ended early"
                )
            if expected_line != observed_line:
                raise ValueError(
                    f"{output_path.name}:{line_number}: content differs from "
                    "the deterministic generator"
                )


def copy_or_verify_metadata(
    source_path: Path,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    if not source_path.exists():
        raise FileNotFoundError(f"Missing source metadata: {source_path}")
    if args.verify_existing:
        if not output_path.exists():
            raise FileNotFoundError(f"Missing generated metadata: {output_path}")
        if sha256_file(source_path) != sha256_file(output_path):
            raise ValueError(
                f"{output_path.name}: metadata differs from tracked source"
            )
        return

    if output_path.exists() and not args.overwrite:
        raise FileExistsError(
            f"{output_path} already exists; use --verify-existing to check it "
            "or --overwrite to replace it"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".copying.tmp")
    try:
        shutil.copyfile(source_path, temporary_path)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def write_manifest(
    records: list[dict[str, object]],
    manifest_path: Path,
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = manifest_path.with_suffix(
        manifest_path.suffix + ".writing.tmp"
    )
    fieldnames = [
        "file",
        "grid_size",
        "cells",
        "source_genes",
        "generated_genes",
        "copies_per_gene",
        "seed",
        "relative_noise_sd",
        "minimum_noise_sd",
        "source_sha256",
        "generated_sha256",
        "generated_bytes",
        "metadata_sha256",
        "numpy_version",
    ]
    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(records)
        os.replace(temporary_path, manifest_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def main() -> None:
    args = parse_args()
    validate_args(args)
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    manifest_path = args.manifest.resolve()
    if source_dir == output_dir:
        raise ValueError(
            "Source and output directories must differ so tracked inputs cannot "
            "be overwritten."
        )
    records: list[dict[str, object]] = []

    for source_path in discover_sources(source_dir):
        grid = grid_size(source_path)
        output_path = output_dir / source_path.name
        source_meta = source_dir / f"sccube_meta_{grid}.csv"
        output_meta = output_dir / source_meta.name

        if args.verify_existing:
            verify_data_file(source_path, output_path, args)
            action = "verified"
        else:
            generate_data_file(source_path, output_path, args)
            action = "generated"
        copy_or_verify_metadata(source_meta, output_meta, args)

        records.append(
            {
                "file": output_path.name,
                "grid_size": grid,
                "cells": grid * grid,
                "source_genes": args.source_genes,
                "generated_genes": args.source_genes * args.copies_per_gene,
                "copies_per_gene": args.copies_per_gene,
                "seed": args.base_seed + grid,
                "relative_noise_sd": args.relative_noise_sd,
                "minimum_noise_sd": args.minimum_noise_sd,
                "source_sha256": sha256_file(source_path),
                "generated_sha256": sha256_file(output_path),
                "generated_bytes": output_path.stat().st_size,
                "metadata_sha256": sha256_file(output_meta),
                "numpy_version": np.__version__,
            }
        )
        print(
            f"{action} {output_path.name}: "
            f"{args.source_genes * args.copies_per_gene} genes x "
            f"{grid * grid} cells"
        )

    write_manifest(records, manifest_path)
    print(f"wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
