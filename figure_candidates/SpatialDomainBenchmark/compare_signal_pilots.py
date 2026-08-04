#!/usr/bin/env python3
"""Compare paired pilot runs across domain-expression fold changes."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from evaluate_pilot import DISPLAY, METHODS, evaluate


SCRIPT_DIR = Path(__file__).resolve().parent
DATASETS = (
    "irregular_mrf_seed0_hard",
    "irregular_mrf_seed0_intermediate",
    "irregular_mrf_seed0_moderate",
)
PAIRED_FILES = ("coordinates.tsv", "truth.tsv", "genes.tsv", "gene_truth.tsv")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_pairing(data_root: Path) -> None:
    reference = data_root / DATASETS[0]
    for dataset in DATASETS[1:]:
        candidate = data_root / dataset
        for filename in PAIRED_FILES:
            if sha256(candidate / filename) != sha256(reference / filename):
                raise ValueError(f"Unpaired pilot input: {dataset}/{filename}")
    count_hashes = {sha256(data_root / dataset / "counts.mtx.gz") for dataset in DATASETS}
    if len(count_hashes) != len(DATASETS):
        raise ValueError("Signal pilots must have distinct count matrices")


def metric_table(summary: pd.DataFrame, metric: str) -> str:
    fold_changes = sorted(summary["domain_fold_change"].unique())
    headers = ["Method", *[f"FC {value:.1f}" for value in fold_changes]]
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for method in METHODS:
        subset = summary[summary["method"] == method].set_index("domain_fold_change")
        values = [DISPLAY[method], *[f"{subset.loc[value, metric]:.3f}" for value in fold_changes]]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def render_comparison(summary: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=180, sharex=True, sharey=True)
    colors = plt.get_cmap("tab10").colors
    metrics = (
        ("adjusted_rand_index", "Adjusted Rand index"),
        ("normalized_mutual_information", "Normalized mutual information"),
    )
    for axis, (metric, title) in zip(axes, metrics):
        for index, method in enumerate(METHODS):
            subset = summary[summary["method"] == method].sort_values(
                "domain_fold_change"
            )
            is_baseline = method in {"expression_pca_kmeans", "coordinates_spectral"}
            axis.plot(
                subset["domain_fold_change"],
                subset[metric],
                marker="o",
                linewidth=1.8,
                markersize=4,
                linestyle="--" if is_baseline else "-",
                color=colors[index],
                label=DISPLAY[method],
            )
        axis.set_title(title)
        axis.set_xlabel("Domain-expression fold change")
        axis.set_xticks([1.6, 2.0, 2.5])
        axis.set_ylim(-0.05, 1.05)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    axes[0].set_ylabel("Agreement with held-out domain truth")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Paired spatial-domain pilot signal sweep (layout seed 0)")
    fig.tight_layout(rect=(0, 0.15, 1, 0.94))
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, output: Path) -> None:
    ari = summary.pivot(index="method", columns="domain_fold_change", values="adjusted_rand_index")
    lines = [
        "# Paired pilot signal comparison",
        "",
        "All three pilots use byte-identical coordinates, truth labels, gene identifiers, "
        "and gene annotations. Only the domain-expression fold change and resulting count "
        "matrix differ. This is a single-layout engineering comparison, not the production "
        "benchmark.",
        "",
        "## Engineering observations",
        "",
        f"- The expression baseline remains null at FC 1.6 (ARI {ari.loc['expression_pca_kmeans', 1.6]:.3f}) "
        f"and FC 2.0 ({ari.loc['expression_pca_kmeans', 2.0]:.3f}), then reaches "
        f"ARI {ari.loc['expression_pca_kmeans', 2.5]:.3f} at FC 2.5.",
        f"- SpaGCN changes earlier, from ARI {ari.loc['spagcn', 1.6]:.3f} at FC 1.6 "
        f"to {ari.loc['spagcn', 2.0]:.3f} at FC 2.0.",
        f"- stLearn remains null at FC 2.0 (ARI {ari.loc['stlearn', 2.0]:.3f}) and "
        f"reaches {ari.loc['stlearn', 2.5]:.3f} at FC 2.5.",
        "- GraphST, STAGATE, and spCLUE remain strong and comparatively stable across "
        "all three paired count matrices; BANKSY remains weak.",
        "- These apparent thresholds come from one layout seed and require replication "
        "before biological or method-level interpretation.",
        "",
        "## Adjusted Rand index",
        "",
        metric_table(summary, "adjusted_rand_index"),
        "",
        "## Normalized mutual information",
        "",
        metric_table(summary, "normalized_mutual_information"),
        "",
    ]
    output.write_text("\n".join(lines))


def main() -> None:
    data_root = SCRIPT_DIR / "pilot_data"
    results_root = SCRIPT_DIR / "pilot_results"
    validate_pairing(data_root)
    summary = evaluate(
        data_root,
        results_root,
        ("irregular_mrf",),
        dataset_names=DATASETS,
    )
    if not summary["status"].eq("success").all():
        failed = summary.loc[~summary["status"].eq("success"), ["dataset", "method"]]
        raise RuntimeError(f"Incomplete signal comparison:\n{failed.to_string(index=False)}")
    order = {method: index for index, method in enumerate(METHODS)}
    summary = summary.assign(method_order=summary["method"].map(order)).sort_values(
        ["domain_fold_change", "method_order"]
    )
    summary.drop(columns="method_order").to_csv(
        results_root / "pilot_signal_comparison.tsv", sep="\t", index=False
    )
    render_comparison(summary, results_root / "pilot_signal_comparison.png")
    write_report(summary, results_root / "PILOT_SIGNAL_COMPARISON.md")
    print(
        summary[
            [
                "dataset",
                "domain_fold_change",
                "method_display",
                "adjusted_rand_index",
                "normalized_mutual_information",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
