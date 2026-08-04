#!/usr/bin/env python3
"""Run six named methods and the expression-only production control."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE = SCRIPT_DIR.parents[3]
DEFAULT_DATA_ROOT = SCRIPT_DIR / "experiment_data"
DEFAULT_RESULTS_ROOT = SCRIPT_DIR / "experiment_results"
TOPOLOGIES = ("curved_layers", "irregular_mrf")
DIFFICULTIES = ("hard", "moderate", "easy")
FOLD_CHANGES = {"hard": 2.0, "moderate": 3.0, "easy": 4.0}
LAYOUT_SEEDS = (0, 1, 2)
LEGACY_METHODS = ("graphst", "stagate", "spagcn", "spclue")
NAMED_METHODS = ("graphst", "stagate", "banksy", "spagcn", "spclue", "stlearn")
NEGATIVE_CONTROLS = ("expression_pca_kmeans",)
METHODS = NAMED_METHODS + NEGATIVE_CONTROLS
TOPOLOGY_SEED_OFFSETS = {"curved_layers": 0, "irregular_mrf": 100}


def dataset_name(topology: str, layout_seed: int, difficulty: str) -> str:
    return f"{topology}_seed{layout_seed}_{difficulty}"


def method_seed(topology: str, layout_seed: int) -> int:
    """Pair stochastic method initialization across difficulty levels."""
    return 100000 + TOPOLOGY_SEED_OFFSETS[topology] + layout_seed


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.is_file() else {}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def execute(command: list[str], log_dir: Path, timeout_seconds: float | None) -> int:
    log_dir.mkdir(parents=True, exist_ok=True)
    run_record_path = log_dir / "run.json"
    previous_run_checksum = (
        sha256_file(run_record_path) if run_record_path.is_file() else None
    )
    (log_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    try:
        with (log_dir / "stdout.log").open("w") as stdout, (
            log_dir / "stderr.log"
        ).open("w") as stderr:
            result = subprocess.run(
                command,
                cwd=WORKSPACE,
                stdout=stdout,
                stderr=stderr,
                text=True,
                check=False,
                timeout=timeout_seconds,
            )
        returncode = result.returncode
        failure_status = "orchestrator_failed"
    except subprocess.TimeoutExpired as error:
        returncode = 124
        failure_status = "timed_out"
        with (log_dir / "stderr.log").open("a") as stderr:
            stderr.write(
                f"\nOrchestrator timeout after {float(error.timeout):.1f} seconds.\n"
            )
    current_run_checksum = (
        sha256_file(run_record_path) if run_record_path.is_file() else None
    )
    if returncode != 0 and current_run_checksum == previous_run_checksum:
        (log_dir / "orchestrator_failure.json").write_text(
            json.dumps(
                {
                    "status": failure_status,
                    "returncode": returncode,
                    "command": command,
                    "finished_utc": datetime.now(timezone.utc).isoformat(),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    return returncode


def validate_dataset(
    dataset_dir: Path,
    topology: str,
    difficulty: str,
    layout_seed: int,
) -> None:
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing dataset manifest: {manifest_path}")
    manifest = read_json(manifest_path)
    expected = {
        "dataset": dataset_dir.name,
        "topology": topology,
        "difficulty": difficulty,
        "layout_seed": layout_seed,
        "domain_fold_change": FOLD_CHANGES[difficulty],
        "n_locations": 3000,
        "n_genes": 2000,
        "n_domains": 4,
        "n_cell_types": 10,
    }
    for field, value in expected.items():
        if manifest.get(field) != value:
            raise ValueError(
                f"{manifest_path}: expected {field}={value!r}, "
                f"found {manifest.get(field)!r}"
            )
    for filename in (
        "counts.mtx.gz",
        "latent_mean.npy.gz",
        "genes.tsv",
        "coordinates.tsv",
        "truth.tsv",
        "gene_truth.tsv",
    ):
        if not (dataset_dir / filename).is_file():
            raise FileNotFoundError(f"Missing dataset file: {dataset_dir / filename}")


def successful_output(run_dir: Path, dataset: str, method: str) -> bool:
    record = read_json(run_dir / "run.json")
    assignments = run_dir / "assignments.tsv"
    if not (
        record.get("status") == "success"
        and record.get("dataset") == dataset
        and record.get("method") == method
        and assignments.is_file()
    ):
        return False
    with assignments.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 3000:
        return False
    location_ids = [row.get("location_id", "") for row in rows]
    predicted = [row.get("predicted_domain", "") for row in rows]
    return (
        len(set(location_ids)) == 3000
        and "" not in location_ids
        and "" not in predicted
        and len(set(predicted)) == 4
    )


def summary_has_failures(path: Path) -> bool:
    if not path.is_file():
        return True
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return not rows or any(row.get("status") != "success" for row in rows)


def method_command(
    method: str,
    dataset_dir: Path,
    output_dir: Path,
    seed: int,
) -> list[str]:
    if method in LEGACY_METHODS:
        return [
            "conda",
            "run",
            "--no-capture-output",
            "-n",
            "spatial-domain-benchmark",
            "python",
            str(SCRIPT_DIR / "run_legacy_method.py"),
            "--method",
            method,
            "--dataset",
            str(dataset_dir),
            "--output",
            str(output_dir),
            "--seed",
            str(seed),
        ]
    if method == "stlearn":
        return [
            "conda",
            "run",
            "--no-capture-output",
            "-n",
            "spatial-domain-benchmark-stlearn",
            "python",
            str(SCRIPT_DIR / "run_stlearn.py"),
            "--dataset",
            str(dataset_dir),
            "--output",
            str(output_dir),
            "--seed",
            str(seed),
        ]
    if method == "banksy":
        return [
            "Rscript",
            str(SCRIPT_DIR / "run_banksy.R"),
            str(dataset_dir),
            str(output_dir),
            str(seed),
        ]
    if method == "expression_pca_kmeans":
        return [
            "conda",
            "run",
            "--no-capture-output",
            "-n",
            "simspace-repro",
            "python",
            str(SCRIPT_DIR / "run_baselines.py"),
            "--dataset",
            str(dataset_dir),
            "--output-root",
            str(output_dir.parent),
            "--seed",
            str(seed),
            "--baselines",
            method,
        ]
    raise ValueError(f"Unknown method: {method}")


def evaluation_command(
    data_root: Path,
    results_root: Path,
    topologies: tuple[str, ...],
    difficulties: tuple[str, ...],
    layout_seeds: tuple[int, ...],
    methods: tuple[str, ...],
) -> list[str]:
    return [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "simspace-repro",
        "python",
        str(SCRIPT_DIR / "evaluate_experiment.py"),
        "--data-root",
        str(data_root),
        "--results-root",
        str(results_root),
        "--topologies",
        *topologies,
        "--difficulties",
        *difficulties,
        "--layout-seeds",
        *(str(seed) for seed in layout_seeds),
        "--methods",
        *methods,
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate, run, and aggregate the production benchmark."
    )
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--skip-generation", action="store_true")
    parser.add_argument("--resume-successful", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--timeout-minutes", type=float, default=90.0,
        help="Per-method timeout. Generation and evaluation are not time limited.",
    )
    parser.add_argument(
        "--topologies", nargs="+", choices=TOPOLOGIES, default=TOPOLOGIES
    )
    parser.add_argument(
        "--difficulties", nargs="+", choices=DIFFICULTIES, default=DIFFICULTIES
    )
    parser.add_argument(
        "--layout-seeds", nargs="+", type=int, choices=LAYOUT_SEEDS, default=LAYOUT_SEEDS
    )
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    args = parser.parse_args()
    if args.timeout_minutes <= 0:
        raise ValueError("--timeout-minutes must be positive")
    if args.resume_successful and not args.skip_generation:
        raise ValueError(
            "--resume-successful requires --skip-generation so completed runs "
            "cannot be paired with regenerated inputs"
        )

    topologies = tuple(dict.fromkeys(args.topologies))
    difficulties = tuple(dict.fromkeys(args.difficulties))
    layout_seeds = tuple(dict.fromkeys(args.layout_seeds))
    methods = tuple(dict.fromkeys(args.methods))
    data_root = args.data_root.resolve()
    results_root = args.results_root.resolve()
    dataset_specs = [
        (topology, difficulty, seed)
        for topology in topologies
        for seed in layout_seeds
        for difficulty in difficulties
    ]
    expected_runs = len(dataset_specs) * len(methods)

    generation = [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "simspace-repro",
        "python",
        str(SCRIPT_DIR / "generate_experiment.py"),
        "--output-root",
        str(data_root),
        "--topologies",
        *topologies,
        "--difficulties",
        *difficulties,
        "--layout-seeds",
        *(str(seed) for seed in layout_seeds),
    ]
    run_specs: list[tuple[str, str, Path, list[str]]] = []
    for topology, difficulty, seed in dataset_specs:
        name = dataset_name(topology, seed, difficulty)
        dataset_dir = data_root / name
        for method in methods:
            output_dir = results_root / name / method
            run_specs.append(
                (
                    name,
                    method,
                    output_dir,
                    method_command(method, dataset_dir, output_dir, method_seed(topology, seed)),
                )
            )
    evaluate = evaluation_command(
        data_root, results_root, topologies, difficulties, layout_seeds, methods
    )

    print(
        f"PLAN {len(dataset_specs)} datasets x {len(methods)} methods = "
        f"{expected_runs} runs",
        flush=True,
    )
    if args.dry_run:
        if not args.skip_generation:
            print("GENERATE " + shlex.join(generation))
        for name, method, _, command in run_specs:
            print(f"RUN {name} / {method}: {shlex.join(command)}")
        print("EVALUATE " + shlex.join(evaluate))
        return

    existing_successes = [
        f"{name} / {method}"
        for name, method, output_dir, _ in run_specs
        if successful_output(output_dir, name, method)
    ]
    if existing_successes and not args.resume_successful:
        raise FileExistsError(
            "Validated successful outputs already exist. Use a new results root, "
            "or use --skip-generation --resume-successful. First existing result: "
            + existing_successes[0]
        )

    results_root.mkdir(parents=True, exist_ok=True)
    plan = {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "topologies": list(topologies),
        "difficulties": {difficulty: FOLD_CHANGES[difficulty] for difficulty in difficulties},
        "layout_seeds": list(layout_seeds),
        "methods": list(methods),
        "named_methods": [method for method in methods if method in NAMED_METHODS],
        "negative_controls": [
            method for method in methods if method in NEGATIVE_CONTROLS
        ],
        "n_datasets": len(dataset_specs),
        "n_planned_runs": expected_runs,
        "method_seed_policy": "100000 + fixed topology offset + layout seed; paired across difficulty",
        "method_seeds": {
            f"{topology}_seed{seed}": method_seed(topology, seed)
            for topology in topologies
            for seed in layout_seeds
        },
        "timeout_minutes": args.timeout_minutes,
        "resume_successful": args.resume_successful,
        "implementation_checksums": {
            filename: sha256_file(SCRIPT_DIR / filename)
            for filename in (
                "benchmark-config.draft.yml",
                "environment.yml",
                "environment-stlearn.yml",
                "generate_pilot.py",
                "generate_experiment.py",
                "pilot_common.py",
                "run_legacy_method.py",
                "run_banksy.R",
                "run_baselines.py",
                "run_stlearn.py",
                "run_experiment.py",
                "evaluate_experiment.py",
            )
        },
    }
    (results_root / "experiment_plan.json").write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n"
    )

    if not args.skip_generation:
        print("GENERATE production datasets", flush=True)
        code = execute(generation, results_root / "generation", timeout_seconds=None)
        if code != 0:
            raise SystemExit(code)

    for topology, difficulty, seed in dataset_specs:
        validate_dataset(
            data_root / dataset_name(topology, seed, difficulty),
            topology,
            difficulty,
            seed,
        )

    run_failures: list[str] = []
    skipped_successes = 0
    timeout_seconds = 60.0 * args.timeout_minutes
    for name, method, output_dir, command in run_specs:
        if args.resume_successful and successful_output(output_dir, name, method):
            print(f"SKIP successful {name} / {method}", flush=True)
            skipped_successes += 1
            continue
        print(f"RUN {name} / {method}", flush=True)
        code = execute(command, output_dir, timeout_seconds=timeout_seconds)
        if code != 0 or not successful_output(output_dir, name, method):
            run_failures.append(f"{name} / {method}")

    print(
        f"EVALUATE {expected_runs} expected rows "
        f"({skipped_successes} successful runs resumed)",
        flush=True,
    )
    evaluation_code = execute(
        evaluate, results_root / "evaluation", timeout_seconds=None
    )
    if evaluation_code != 0:
        raise SystemExit(evaluation_code)
    invalid_summary = summary_has_failures(results_root / "experiment_summary.tsv")
    if run_failures or invalid_summary:
        if run_failures:
            print("FAILED RUNS: " + ", ".join(run_failures), file=sys.stderr)
        if invalid_summary:
            print(
                "The evaluation summary contains missing, failed, or invalid outputs.",
                file=sys.stderr,
            )
        raise SystemExit(1)


if __name__ == "__main__":
    main()
