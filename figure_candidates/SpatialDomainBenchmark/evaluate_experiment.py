#!/usr/bin/env python3
"""Aggregate held-out metrics for the production spatial-domain benchmark."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from pilot_common import boundary_f1, boundary_fraction, sha256_file


SCRIPT_DIR = Path(__file__).resolve().parent
TOPOLOGIES = ("curved_layers", "irregular_mrf")
DIFFICULTIES = ("hard", "moderate", "easy")
FOLD_CHANGES = {"hard": 2.0, "moderate": 3.0, "easy": 4.0}
LAYOUT_SEEDS = (0, 1, 2)
NAMED_METHODS = ("graphst", "stagate", "banksy", "spagcn", "spclue", "stlearn")
NEGATIVE_CONTROLS = ("expression_pca_kmeans",)
METHODS = NAMED_METHODS + NEGATIVE_CONTROLS
DISPLAY = {
    "graphst": "GraphST",
    "stagate": "STAGATE",
    "banksy": "BANKSY",
    "spagcn": "SpaGCN",
    "spclue": "spCLUE",
    "stlearn": "stLearn",
    "expression_pca_kmeans": "Expression-only control",
}
METRIC_COLUMNS = (
    "adjusted_rand_index",
    "normalized_mutual_information",
    "boundary_precision",
    "boundary_recall",
    "boundary_f1",
    "produced_cluster_count",
)


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.is_file() else {}


def dataset_name(topology: str, layout_seed: int, difficulty: str) -> str:
    return f"{topology}_seed{layout_seed}_{difficulty}"


def align_assignments(
    truth: pd.DataFrame, assignments: pd.DataFrame
) -> pd.DataFrame:
    required = {"location_id", "predicted_domain"}
    missing_columns = required - set(assignments.columns)
    if missing_columns:
        raise ValueError(
            "Assignments lack required columns: " + ", ".join(sorted(missing_columns))
        )
    if assignments["location_id"].duplicated().any():
        raise ValueError("Assignments contain duplicate location identifiers")
    aligned = truth.merge(
        assignments[["location_id", "predicted_domain"]],
        on="location_id",
        how="left",
        validate="one_to_one",
    )
    if aligned["predicted_domain"].isna().any():
        missing = int(aligned["predicted_domain"].isna().sum())
        raise ValueError(f"Assignments are missing {missing} locations")
    extras = set(assignments["location_id"]) - set(truth["location_id"])
    if extras:
        raise ValueError(f"Assignments contain {len(extras)} unknown locations")
    return aligned


def validate_manifest(
    dataset_dir: Path,
    topology: str,
    difficulty: str,
    layout_seed: int,
) -> dict[str, object]:
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing manifest: {manifest_path}")
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
    for filename, expected_checksum in manifest.get("checksums", {}).items():
        path = dataset_dir / filename
        if not path.is_file() or sha256_file(path) != expected_checksum:
            raise ValueError(f"Missing or checksum-invalid dataset file: {path}")
    return manifest


def initial_status(run_dir: Path) -> tuple[dict[str, object], str]:
    run_path = run_dir / "run.json"
    failure_path = run_dir / "orchestrator_failure.json"
    record = read_json(run_path)
    orchestrator_failure = read_json(failure_path)
    failure_is_newer = bool(orchestrator_failure) and (
        not run_path.is_file()
        or failure_path.stat().st_mtime_ns > run_path.stat().st_mtime_ns
    )
    if failure_is_newer:
        return orchestrator_failure, str(
            orchestrator_failure.get("status", "orchestrator_failed")
        )
    if record:
        return record, str(record.get("status", "unknown"))
    if orchestrator_failure:
        return orchestrator_failure, str(
            orchestrator_failure.get("status", "orchestrator_failed")
        )
    return {}, "missing"


def evaluate(
    data_root: Path,
    results_root: Path,
    topologies: tuple[str, ...],
    difficulties: tuple[str, ...],
    layout_seeds: tuple[int, ...],
    methods: tuple[str, ...],
) -> pd.DataFrame:
    root_manifest_path = data_root / "experiment_manifest.json"
    if not root_manifest_path.is_file():
        raise FileNotFoundError(f"Missing experiment manifest: {root_manifest_path}")
    root_manifest = read_json(root_manifest_path)
    generated_names = {
        str(record["dataset"]) for record in root_manifest.get("datasets", [])
    }
    rows: list[dict[str, object]] = []
    for topology in topologies:
        for layout_seed in layout_seeds:
            for difficulty in difficulties:
                name = dataset_name(topology, layout_seed, difficulty)
                if name not in generated_names:
                    raise ValueError(f"Dataset absent from experiment manifest: {name}")
                dataset_dir = data_root / name
                manifest = validate_manifest(
                    dataset_dir, topology, difficulty, layout_seed
                )
                truth = pd.read_csv(dataset_dir / "truth.tsv", sep="\t")
                coordinates = pd.read_csv(dataset_dir / "coordinates.tsv", sep="\t")
                if not truth["location_id"].equals(coordinates["location_id"]):
                    raise ValueError(f"Truth/coordinate location order differs for {name}")
                xy = coordinates[["x", "y"]].to_numpy(dtype=float)
                true_labels = truth["domain_truth"].to_numpy()
                diagnostics = {
                    "truth_boundary_fraction": boundary_fraction(xy, true_labels),
                    "domain_cell_type_nmi": normalized_mutual_info_score(
                        true_labels, truth["cell_type_truth"]
                    ),
                }
                for method in methods:
                    run_dir = results_root / name / method
                    record, status = initial_status(run_dir)
                    row: dict[str, object] = {
                        "dataset": name,
                        "topology": topology,
                        "difficulty": difficulty,
                        "domain_fold_change": manifest["domain_fold_change"],
                        "layout_seed": layout_seed,
                        "method": method,
                        "method_display": DISPLAY[method],
                        "method_seed": record.get("method_seed", np.nan),
                        "status": status,
                        "wall_time_seconds": record.get("wall_time_seconds", np.nan),
                        "peak_rss_mb": record.get("peak_rss_mb", np.nan),
                        "error_type": record.get("error_type", ""),
                        "error_message": record.get("error_message", ""),
                        **diagnostics,
                    }
                    row.update({column: np.nan for column in METRIC_COLUMNS})
                    assignments_path = run_dir / "assignments.tsv"
                    if status == "success":
                        try:
                            if not assignments_path.is_file():
                                raise FileNotFoundError("assignments.tsv is missing")
                            assignments = pd.read_csv(
                                assignments_path, sep="\t", dtype=str
                            )
                            aligned = align_assignments(truth, assignments)
                            predicted = aligned["predicted_domain"].to_numpy()
                            cluster_count = int(pd.Series(predicted).nunique())
                            if cluster_count != 4:
                                raise ValueError(
                                    f"Expected four nonempty clusters, found {cluster_count}"
                                )
                            precision, recall, f1 = boundary_f1(
                                xy, true_labels, predicted
                            )
                            row.update(
                                {
                                    "adjusted_rand_index": adjusted_rand_score(
                                        true_labels, predicted
                                    ),
                                    "normalized_mutual_information": normalized_mutual_info_score(
                                        true_labels, predicted
                                    ),
                                    "boundary_precision": precision,
                                    "boundary_recall": recall,
                                    "boundary_f1": f1,
                                    "produced_cluster_count": cluster_count,
                                }
                            )
                        except Exception as error:
                            row["status"] = "invalid_output"
                            row["error_type"] = type(error).__name__
                            row["error_message"] = str(error)
                    rows.append(row)
    summary = pd.DataFrame(rows)
    expected_rows = (
        len(topologies) * len(difficulties) * len(layout_seeds) * len(methods)
    )
    if len(summary) != expected_rows:
        raise RuntimeError(f"Expected {expected_rows} summary rows, found {len(summary)}")
    if summary[["dataset", "method"]].duplicated().any():
        raise RuntimeError("Duplicate dataset/method rows in production summary")
    return summary


def format_metric(value: float) -> str:
    return "" if pd.isna(value) else f"{float(value):.3f}"


def write_report(
    summary: pd.DataFrame,
    output: Path,
    topologies: tuple[str, ...],
    difficulties: tuple[str, ...],
    layout_seeds: tuple[int, ...],
    methods: tuple[str, ...],
) -> None:
    named = summary[summary["method"].isin(NAMED_METHODS)]
    controls = summary[summary["method"].isin(NEGATIVE_CONTROLS)]
    successful_named = named[named["status"] == "success"]
    successful_controls = controls[controls["status"] == "success"]
    failures = summary[summary["status"] != "success"]
    lines = [
        "# Spatial-domain production benchmark report",
        "",
        f"Named-method completion: **{len(successful_named)}/{len(named)} runs** across "
        f"{summary['dataset'].nunique()} datasets.",
        "",
    ]
    if len(controls):
        lines.extend(
            [
                f"Expression-only negative-control completion: "
                f"**{len(successful_controls)}/{len(controls)} runs**.",
                "",
            ]
        )
    lines.extend([
        "Each row below summarizes independent layout seeds. Difficulty levels are "
        "paired within pattern and seed; only the nominal domain-expression fold "
        "change differs within each trio.",
        "",
        "## Accuracy and runtime summary",
        "",
        "| Pattern | Difficulty | Method | Completed | Median ARI | Median NMI | Median boundary F1 | Median seconds |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ])
    for topology in topologies:
        for difficulty in difficulties:
            for method in methods:
                group = summary[
                    (summary["topology"] == topology)
                    & (summary["difficulty"] == difficulty)
                    & (summary["method"] == method)
                ]
                complete = group[group["status"] == "success"]
                medians = complete[
                    [
                        "adjusted_rand_index",
                        "normalized_mutual_information",
                        "boundary_f1",
                        "wall_time_seconds",
                    ]
                ].median(numeric_only=True)
                lines.append(
                    "| "
                    + " | ".join(
                        [
                            topology.replace("_", " "),
                            difficulty,
                            DISPLAY[method],
                            f"{len(complete)}/{len(layout_seeds)}",
                            format_metric(medians.get("adjusted_rand_index", np.nan)),
                            format_metric(
                                medians.get("normalized_mutual_information", np.nan)
                            ),
                            format_metric(medians.get("boundary_f1", np.nan)),
                            (
                                ""
                                if pd.isna(medians.get("wall_time_seconds", np.nan))
                                else f"{float(medians['wall_time_seconds']):.1f}"
                            ),
                        ]
                    )
                    + " |"
                )
    lines.extend(["", "## Failures and invalid outputs", ""])
    if failures.empty:
        lines.append("No failures or invalid outputs.")
    else:
        lines.extend(
            [
                "| Dataset | Method | Status | Error |",
                "|---|---|---|---|",
            ]
        )
        for _, row in failures.iterrows():
            error = f"{row.get('error_type', '')}: {row.get('error_message', '')}".strip(
                ": "
            )
            error = error.replace("|", "\\|").replace("\n", " ")
            lines.append(
                f"| {row['dataset']} | {row['method_display']} | "
                f"{row['status']} | {error} |"
            )
    lines.extend(
        [
            "",
            "## Evaluation contract",
            "",
            "- ARI and NMI compare predictions with held-out SimSpace domain labels.",
            "- Boundary F1 uses the fixed six-nearest-neighbor, one-hop-tolerant definition.",
            "- A completed result must contain exactly 3,000 aligned assignments and four nonempty clusters.",
            "- Accuracy metrics are omitted for failures and invalid outputs; completion remains explicit.",
            "- `truth.tsv` is read only during this aggregation step, never by a method runner.",
            "- The expression-only control receives counts but no spatial coordinates in its computation and is reported separately from the six named methods.",
            "",
        ]
    )
    output.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate selected cells of the frozen production design."
    )
    parser.add_argument(
        "--data-root", type=Path, default=SCRIPT_DIR / "experiment_data"
    )
    parser.add_argument(
        "--results-root", type=Path, default=SCRIPT_DIR / "experiment_results"
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

    topologies = tuple(dict.fromkeys(args.topologies))
    difficulties = tuple(dict.fromkeys(args.difficulties))
    layout_seeds = tuple(dict.fromkeys(args.layout_seeds))
    methods = tuple(dict.fromkeys(args.methods))
    data_root = args.data_root.resolve()
    results_root = args.results_root.resolve()
    results_root.mkdir(parents=True, exist_ok=True)
    summary = evaluate(
        data_root, results_root, topologies, difficulties, layout_seeds, methods
    )
    summary.to_csv(results_root / "experiment_summary.tsv", sep="\t", index=False)
    write_report(
        summary,
        results_root / "EXPERIMENT_REPORT.md",
        topologies,
        difficulties,
        layout_seeds,
        methods,
    )
    metadata = {
        "schema_version": "1.0",
        "evaluated_utc": datetime.now(timezone.utc).isoformat(),
        "topologies": list(topologies),
        "difficulties": list(difficulties),
        "layout_seeds": list(layout_seeds),
        "methods": list(methods),
        "named_methods": [method for method in methods if method in NAMED_METHODS],
        "negative_controls": [
            method for method in methods if method in NEGATIVE_CONTROLS
        ],
        "n_datasets": int(summary["dataset"].nunique()),
        "n_expected_runs": len(summary),
        "n_successful_runs": int((summary["status"] == "success").sum()),
        "summary_sha256": sha256_file(results_root / "experiment_summary.tsv"),
    }
    (results_root / "evaluation_manifest.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Wrote {len(summary)} rows for {summary['dataset'].nunique()} datasets; "
        f"{metadata['n_successful_runs']} runs completed successfully."
    )


if __name__ == "__main__":
    main()
