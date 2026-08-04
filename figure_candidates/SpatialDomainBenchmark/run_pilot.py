#!/usr/bin/env python3
"""Orchestrate the single-MRF pilot across its three runtimes."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_ROOT = SCRIPT_DIR / "pilot_data"
RESULTS_ROOT = SCRIPT_DIR / "pilot_results"
LEGACY_METHODS = ("graphst", "stagate", "spagcn", "spclue")
ALL_METHODS = LEGACY_METHODS + ("banksy", "stlearn", "baselines")
PILOT_DIFFICULTY = "intermediate"
PILOT_DOMAIN_FOLD_CHANGE = 3.0


def execute(command: list[str], log_dir: Path) -> int:
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    with (log_dir / "stdout.log").open("w") as stdout, (
        log_dir / "stderr.log"
    ).open("w") as stderr:
        result = subprocess.run(
            command,
            cwd=SCRIPT_DIR.parents[3],
            stdout=stdout,
            stderr=stderr,
            text=True,
            check=False,
        )
    if result.returncode != 0 and not (log_dir / "run.json").exists():
        (log_dir / "orchestrator_failure.json").write_text(
            json.dumps({"returncode": result.returncode, "command": command}, indent=2)
            + "\n"
        )
    return result.returncode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-generation", action="store_true")
    parser.add_argument("--methods", nargs="+", choices=ALL_METHODS, default=ALL_METHODS)
    parser.add_argument("--difficulty-label", default=PILOT_DIFFICULTY)
    parser.add_argument(
        "--domain-fold-change", type=float, default=PILOT_DOMAIN_FOLD_CHANGE
    )
    args = parser.parse_args()
    if not args.difficulty_label.replace("_", "").isalnum():
        raise ValueError("difficulty label must contain only letters, digits, or underscores")
    pilot_datasets = (f"irregular_mrf_seed0_{args.difficulty_label}",)
    if not args.skip_generation:
        code = execute(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "simspace-repro",
                "python",
                str(SCRIPT_DIR / "generate_pilot.py"),
                "--topologies",
                "irregular_mrf",
                "--difficulty-label",
                args.difficulty_label,
                "--domain-fold-change",
                str(args.domain_fold_change),
            ],
            RESULTS_ROOT / "generation",
        )
        if code != 0:
            raise SystemExit(code)

    datasets = [DATA_ROOT / name for name in pilot_datasets]
    missing = [str(path) for path in datasets if not path.is_dir()]
    if missing:
        raise FileNotFoundError("Missing pilot dataset(s): " + ", ".join(missing))
    for dataset in datasets:
        manifest = json.loads((dataset / "manifest.json").read_text())
        if manifest.get("difficulty") != args.difficulty_label:
            raise ValueError(f"Difficulty mismatch in {dataset / 'manifest.json'}")
        if float(manifest.get("domain_fold_change")) != args.domain_fold_change:
            raise ValueError(f"Domain-fold-change mismatch in {dataset / 'manifest.json'}")
    run_failures: list[str] = []
    for dataset_index, dataset in enumerate(datasets):
        seed = 100000 + dataset_index
        dataset_result = RESULTS_ROOT / dataset.name
        for method in args.methods:
            if method in LEGACY_METHODS:
                output = dataset_result / method
                command = [
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
                    str(dataset),
                    "--output",
                    str(output),
                    "--seed",
                    str(seed),
                ]
                print(f"RUN {dataset.name} / {method}", flush=True)
                if execute(command, output) != 0:
                    run_failures.append(f"{dataset.name} / {method}")
            elif method == "stlearn":
                output = dataset_result / method
                command = [
                    "conda",
                    "run",
                    "--no-capture-output",
                    "-n",
                    "spatial-domain-benchmark-stlearn",
                    "python",
                    str(SCRIPT_DIR / "run_stlearn.py"),
                    "--dataset",
                    str(dataset),
                    "--output",
                    str(output),
                    "--seed",
                    str(seed),
                ]
                print(f"RUN {dataset.name} / stlearn", flush=True)
                if execute(command, output) != 0:
                    run_failures.append(f"{dataset.name} / stlearn")
            elif method == "banksy":
                output = dataset_result / method
                command = [
                    "Rscript",
                    str(SCRIPT_DIR / "run_banksy.R"),
                    str(dataset),
                    str(output),
                    str(seed),
                ]
                print(f"RUN {dataset.name} / banksy", flush=True)
                if execute(command, output) != 0:
                    run_failures.append(f"{dataset.name} / banksy")
            elif method == "baselines":
                command = [
                    "conda",
                    "run",
                    "--no-capture-output",
                    "-n",
                    "simspace-repro",
                    "python",
                    str(SCRIPT_DIR / "run_baselines.py"),
                    "--dataset",
                    str(dataset),
                    "--output-root",
                    str(dataset_result),
                    "--seed",
                    str(seed),
                ]
                print(f"RUN {dataset.name} / baselines", flush=True)
                if execute(command, dataset_result / "baselines") != 0:
                    run_failures.append(f"{dataset.name} / baselines")

    evaluate_command = [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        "simspace-repro",
        "python",
        str(SCRIPT_DIR / "evaluate_pilot.py"),
        "--topologies",
        "irregular_mrf",
        "--datasets",
        *pilot_datasets,
    ]
    code = execute(evaluate_command, RESULTS_ROOT / "evaluation")
    if code != 0:
        raise SystemExit(code)
    if run_failures:
        print("FAILED RUNS: " + ", ".join(run_failures), file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
