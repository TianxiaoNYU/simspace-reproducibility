#!/usr/bin/env python3
"""Evaluate and render all completed pilot runs from held-out truth."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from pilot_common import boundary_f1, boundary_fraction


NAMED_METHODS = (
    "graphst",
    "stagate",
    "banksy",
    "spagcn",
    "spclue",
    "stlearn",
)
METHODS = NAMED_METHODS + (
    "expression_pca_kmeans",
    "coordinates_spectral",
)
DISPLAY = {
    "graphst": "GraphST",
    "stagate": "STAGATE",
    "banksy": "BANKSY",
    "spagcn": "SpaGCN",
    "spclue": "spCLUE",
    "stlearn": "stLearn",
    "expression_pca_kmeans": "Expression baseline",
    "coordinates_spectral": "Coordinate baseline",
}


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text()) if path.exists() else {}


def align_assignments(
    truth: pd.DataFrame, assignments: pd.DataFrame
) -> pd.DataFrame:
    if assignments["location_id"].duplicated().any():
        raise ValueError("Duplicate location identifiers in assignments")
    aligned = truth.merge(assignments, on="location_id", how="left", validate="one_to_one")
    if aligned["predicted_domain"].isna().any():
        missing = int(aligned["predicted_domain"].isna().sum())
        raise ValueError(f"Assignments are missing {missing} locations")
    extra = set(assignments["location_id"]) - set(truth["location_id"])
    if extra:
        raise ValueError(f"Assignments contain {len(extra)} unknown locations")
    return aligned


def selected_datasets(
    data_root: Path,
    topologies: tuple[str, ...],
    dataset_names: tuple[str, ...] = (),
) -> list[Path]:
    datasets = []
    candidates = (
        [data_root / name for name in dataset_names]
        if dataset_names
        else sorted(candidate for candidate in data_root.iterdir() if candidate.is_dir())
    )
    for path in candidates:
        if not path.is_dir():
            raise FileNotFoundError(f"Missing dataset directory: {path}")
        if read_json(path / "manifest.json").get("topology") in topologies:
            datasets.append(path)
    if not datasets:
        raise ValueError(f"No datasets found for topologies: {', '.join(topologies)}")
    return datasets


def evaluate(
    data_root: Path,
    results_root: Path,
    topologies: tuple[str, ...],
    dataset_names: tuple[str, ...] = (),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_dir in selected_datasets(data_root, topologies, dataset_names):
        manifest = read_json(dataset_dir / "manifest.json")
        niche_mrf = manifest.get("niche_mrf") or {}
        truth = pd.read_csv(dataset_dir / "truth.tsv", sep="\t")
        coordinates = pd.read_csv(dataset_dir / "coordinates.tsv", sep="\t")
        xy = coordinates[["x", "y"]].to_numpy(dtype=float)
        domain_cell_type_nmi = normalized_mutual_info_score(
            truth["domain_truth"], truth["cell_type_truth"]
        )
        truth_boundary_fraction = boundary_fraction(
            xy, truth["domain_truth"].to_numpy()
        )
        for method in METHODS:
            run_dir = results_root / dataset_dir.name / method
            record = read_json(run_dir / "run.json")
            row: dict[str, object] = {
                "dataset": dataset_dir.name,
                "topology": manifest.get("topology"),
                "difficulty": manifest.get("difficulty"),
                "domain_fold_change": manifest.get("domain_fold_change", np.nan),
                "method": method,
                "method_display": DISPLAY[method],
                "status": record.get("status", "missing"),
                "wall_time_seconds": record.get("wall_time_seconds", np.nan),
                "peak_rss_mb": record.get("peak_rss_mb", np.nan),
                "error_type": record.get("error_type", ""),
                "error_message": record.get("error_message", ""),
                "truth_boundary_fraction": truth_boundary_fraction,
                "domain_cell_type_nmi": domain_cell_type_nmi,
                "grid_rows": manifest.get("grid_shape", [np.nan, np.nan])[0],
                "grid_columns": manifest.get("grid_shape", [np.nan, np.nan])[1],
                "retention_fraction": manifest.get("retention_fraction", np.nan),
                "coordinate_jitter": manifest.get("coordinate_jitter", np.nan),
                "enriched_probability_mass": manifest.get(
                    "enriched_probability_mass", np.nan
                ),
                "shared_probability_mass": manifest.get(
                    "shared_probability_mass", np.nan
                ),
                "niche_mrf_phi": niche_mrf.get("phi", np.nan),
                "niche_mrf_sweeps": niche_mrf.get("sweeps", np.nan),
                "niche_mrf_neighborhood_radius": niche_mrf.get(
                    "neighborhood_radius", np.nan
                ),
            }
            assignments_path = run_dir / "assignments.tsv"
            if row["status"] == "success" and assignments_path.exists():
                assignments = pd.read_csv(assignments_path, sep="\t", dtype=str)
                aligned = align_assignments(truth, assignments)
                true_labels = aligned["domain_truth"].to_numpy()
                predicted = aligned["predicted_domain"].to_numpy()
                precision, recall, f1 = boundary_f1(xy, true_labels, predicted)
                row.update(
                    {
                        "adjusted_rand_index": adjusted_rand_score(true_labels, predicted),
                        "normalized_mutual_information": normalized_mutual_info_score(
                            true_labels, predicted
                        ),
                        "boundary_precision": precision,
                        "boundary_recall": recall,
                        "boundary_f1": f1,
                        "produced_cluster_count": int(pd.Series(predicted).nunique()),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows)


def render_overview(
    data_root: Path,
    results_root: Path,
    output: Path,
    topologies: tuple[str, ...],
    dataset_names: tuple[str, ...] = (),
) -> None:
    dataset_dirs = selected_datasets(data_root, topologies, dataset_names)
    columns = ("truth",) + METHODS
    fig, axes = plt.subplots(
        len(dataset_dirs), len(columns), figsize=(2.1 * len(columns), 4.2), dpi=180
    )
    palette = np.asarray(
        ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#999999"]
    )
    if len(dataset_dirs) == 1:
        axes = axes[None, :]
    for row_index, dataset_dir in enumerate(dataset_dirs):
        coordinates = pd.read_csv(dataset_dir / "coordinates.tsv", sep="\t")
        truth = pd.read_csv(dataset_dir / "truth.tsv", sep="\t")
        for column_index, method in enumerate(columns):
            ax = axes[row_index, column_index]
            labels: np.ndarray | None
            if method == "truth":
                labels = truth["domain_truth"].to_numpy()
                title = "Truth"
            else:
                path = results_root / dataset_dir.name / method / "assignments.tsv"
                title = DISPLAY[method]
                if path.exists():
                    assignments = pd.read_csv(path, sep="\t", dtype=str)
                    aligned = align_assignments(truth, assignments)
                    labels = pd.factorize(aligned["predicted_domain"], sort=True)[0]
                else:
                    labels = None
            if labels is None:
                ax.text(0.5, 0.5, "FAILED", ha="center", va="center", color="#B2182B")
            else:
                ax.scatter(
                    coordinates["x"],
                    -coordinates["y"],
                    c=palette[np.asarray(labels, dtype=int) % len(palette)],
                    s=1.0,
                    linewidths=0,
                    rasterized=True,
                )
            if row_index == 0:
                ax.set_title(title, fontsize=8)
            if column_index == 0:
                topology = read_json(dataset_dir / "manifest.json")["topology"]
                ax.set_ylabel(topology.replace("_", " "), fontsize=8)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
    difficulties = sorted(
        {str(read_json(path / "manifest.json").get("difficulty")) for path in dataset_dirs}
    )
    fig.suptitle(
        f"Spatial-domain MRF {'/'.join(difficulties)}-signal pilot: "
        "truth and predicted partitions",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def markdown_table(frame: pd.DataFrame) -> str:
    columns = ["Dataset", "Method", "Status", "ARI", "NMI", "Boundary F1", "Seconds"]
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for _, row in frame.iterrows():
        values = [
            str(row["dataset"]),
            str(row["method_display"]),
            str(row["status"]),
            "" if pd.isna(row.get("adjusted_rand_index")) else f"{row['adjusted_rand_index']:.3f}",
            "" if pd.isna(row.get("normalized_mutual_information")) else f"{row['normalized_mutual_information']:.3f}",
            "" if pd.isna(row.get("boundary_f1")) else f"{row['boundary_f1']:.3f}",
            "" if pd.isna(row.get("wall_time_seconds")) else f"{float(row['wall_time_seconds']):.1f}",
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_report(summary: pd.DataFrame, output: Path) -> None:
    named = summary[summary["method"].isin(NAMED_METHODS)]
    successes = int((named["status"] == "success").sum())
    total = len(named)
    failures = named[named["status"] != "success"]
    boundary_diagnostics = summary.drop_duplicates("dataset")[[
        "dataset",
        "difficulty",
        "domain_fold_change",
        "truth_boundary_fraction",
        "domain_cell_type_nmi",
        "grid_rows",
        "grid_columns",
        "retention_fraction",
        "coordinate_jitter",
        "enriched_probability_mass",
        "shared_probability_mass",
        "niche_mrf_phi",
        "niche_mrf_sweeps",
        "niche_mrf_neighborhood_radius",
    ]]
    lines = [
        "# Spatial-domain MRF benchmark pilot report",
        "",
        f"Named-method completion: **{successes}/{total} runs**.",
        "",
        "This is an engineering and design pilot, not the production benchmark. "
        "Metrics are shown to detect broken adapters or trivial simulations; they "
        "must not be used to tune method parameters against truth.",
        "",
        "## Design diagnostics",
        "",
    ]
    for _, row in boundary_diagnostics.iterrows():
        fraction = float(row["truth_boundary_fraction"])
        lines.append(
            f"- {row['dataset']}: {row['difficulty']} signal with domain fold change "
            f"{float(row['domain_fold_change']):.1f}; {int(row['grid_rows'])} x "
            f"{int(row['grid_columns'])} grid, {float(row['retention_fraction']):.1%} "
            f"retention, jitter SD {float(row['coordinate_jitter']):.2f}, and "
            f"{fraction:.1%} truth-boundary nodes in the fixed 6-NN graph; niche "
            f"MRF phi {float(row['niche_mrf_phi']):.2f}, "
            f"{int(row['niche_mrf_sweeps'])} sweeps, Manhattan radius "
            f"{int(row['niche_mrf_neighborhood_radius'])}; cell-type probability "
            f"mass {float(row['enriched_probability_mass']):.0%} enriched / "
            f"{float(row['shared_probability_mass']):.0%} shared, with realized "
            f"domain/cell-type NMI {float(row['domain_cell_type_nmi']):.3f}."
        )
    lines.extend([
        "",
        markdown_table(summary),
        "",
        "## Failures",
        "",
    ])
    if failures.empty:
        lines.append("No named-method failures.")
    else:
        for _, row in failures.iterrows():
            lines.append(
                f"- {row['dataset']} / {row['method_display']}: "
                f"{row.get('error_type', '')} — {row.get('error_message', '')}"
            )
    lines.extend(
        [
            "",
            "## Pilot-specific deviations and checks",
            "",
            "- BANKSY uses its native mclust endpoint with `lambda=0.8` and `G=4`.",
            "- stLearn uses only physical-distance and expression SME weights; the "
            "unused morphology placeholder is checked for output invariance.",
            "- All method metrics are recomputed here from held-out `truth.tsv`; "
            "runner input files contain no truth labels.",
            "- The low-signal pilots show that spatially fragmented predictions can "
            "receive a high one-hop-tolerant boundary F1 despite near-zero ARI/NMI; "
            "boundary F1 must therefore be interpreted jointly with partition metrics.",
            "",
        ]
    )
    output.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    script_dir = Path(__file__).resolve().parent
    parser.add_argument("--data-root", type=Path, default=script_dir / "pilot_data")
    parser.add_argument("--results-root", type=Path, default=script_dir / "pilot_results")
    parser.add_argument(
        "--topologies",
        nargs="+",
        choices=("curved_layers", "irregular_mrf"),
        default=("irregular_mrf",),
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=(),
        help="Optional exact dataset directory names to aggregate.",
    )
    args = parser.parse_args()
    args.results_root.mkdir(parents=True, exist_ok=True)
    topologies = tuple(args.topologies)
    dataset_names = tuple(args.datasets)
    summary = evaluate(args.data_root, args.results_root, topologies, dataset_names)
    summary.to_csv(args.results_root / "pilot_summary.tsv", sep="\t", index=False)
    render_overview(
        args.data_root,
        args.results_root,
        args.results_root / "pilot_overview.png",
        topologies,
        dataset_names,
    )
    write_report(summary, args.results_root / "PILOT_REPORT.md")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
