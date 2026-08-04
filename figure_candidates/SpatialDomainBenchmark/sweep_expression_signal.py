#!/usr/bin/env python3
"""Sweep domain fold change for the frozen expression-only baseline."""

from __future__ import annotations

import gc
import os
import time
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spatial-domain-pilot-matplotlib")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from generate_pilot import N_GENES, N_LOCATIONS, simulate_counts, simulate_locations
from run_baselines import expression_labels


SCRIPT_DIR = Path(__file__).resolve().parent
FOLD_CHANGES = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
LAYOUT_SEEDS = (0, 1, 2)


def run_sweep() -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for layout_seed in LAYOUT_SEEDS:
        _, domain_truth, cell_type_truth = simulate_locations(
            "irregular_mrf", layout_seed
        )
        method_seed = 100000 + layout_seed
        for fold_change in FOLD_CHANGES:
            started = time.perf_counter()
            counts, _ = simulate_counts(
                domain_truth,
                cell_type_truth,
                seed=layout_seed,
                domain_fold_change=fold_change,
            )
            labels, embedding = expression_labels(counts, method_seed)
            rows.append(
                {
                    "layout_seed": layout_seed,
                    "method_seed": method_seed,
                    "domain_fold_change": fold_change,
                    "adjusted_rand_index": adjusted_rand_score(domain_truth, labels),
                    "normalized_mutual_information": normalized_mutual_info_score(
                        domain_truth, labels
                    ),
                    "produced_cluster_count": int(np.unique(labels).size),
                    "n_locations": N_LOCATIONS,
                    "n_genes": N_GENES,
                    "wall_time_seconds": time.perf_counter() - started,
                }
            )
            del counts, labels, embedding
            gc.collect()
            print(
                f"seed={layout_seed} fold_change={fold_change:.1f} "
                f"ARI={rows[-1]['adjusted_rand_index']:.3f}",
                flush=True,
            )
    return pd.DataFrame(rows)


def metric_table(frame: pd.DataFrame, metric: str) -> str:
    pivot = frame.pivot(
        index="domain_fold_change", columns="layout_seed", values=metric
    )
    headers = ["Fold change", *[f"Seed {seed}" for seed in LAYOUT_SEEDS], "Median"]
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for fold_change, row in pivot.iterrows():
        values = [
            f"{fold_change:.1f}",
            *[f"{row[seed]:.3f}" for seed in LAYOUT_SEEDS],
            f"{row.median():.3f}",
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def render_sweep(frame: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=180, sharex=True, sharey=True)
    metrics = (
        ("adjusted_rand_index", "Adjusted Rand index"),
        ("normalized_mutual_information", "Normalized mutual information"),
    )
    seed_colors = ("#56B4E9", "#E69F00", "#009E73")
    for axis, (metric, title) in zip(axes, metrics):
        for seed, color in zip(LAYOUT_SEEDS, seed_colors):
            subset = frame[frame["layout_seed"] == seed].sort_values(
                "domain_fold_change"
            )
            axis.plot(
                subset["domain_fold_change"],
                subset[metric],
                color=color,
                marker="o",
                linewidth=1.2,
                markersize=3.5,
                alpha=0.75,
                label=f"Seed {seed}",
            )
        median = frame.groupby("domain_fold_change")[metric].median()
        axis.plot(
            median.index,
            median.values,
            color="#000000",
            marker="o",
            linewidth=2.2,
            markersize=4,
            label="Median",
        )
        axis.set_title(title)
        axis.set_xlabel("Domain-expression fold change")
        axis.set_xticks(FOLD_CHANGES)
        axis.set_ylim(-0.05, 1.05)
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.6)
    axes[0].set_ylabel("Agreement with held-out domain truth")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Revised heterogeneous-program expression baseline sweep")
    fig.tight_layout(rect=(0, 0.12, 1, 0.94))
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def first_tested_success(frame: pd.DataFrame, seed: int, threshold: float = 0.9) -> str:
    subset = frame[
        (frame["layout_seed"] == seed) & (frame["adjusted_rand_index"] >= threshold)
    ]
    if subset.empty:
        return f">{max(FOLD_CHANGES):.1f}"
    return f"{subset['domain_fold_change'].min():.1f}"


def write_report(frame: pd.DataFrame, output: Path) -> None:
    thresholds = ", ".join(
        f"seed {seed}: {first_tested_success(frame, seed)}" for seed in LAYOUT_SEEDS
    )
    tier_ari = frame.groupby("domain_fold_change")["adjusted_rand_index"].median()
    moderate = frame[frame["domain_fold_change"] == 3.0]["adjusted_rand_index"]
    lines = [
        "# Revised expression-baseline domain-signal sweep",
        "",
        "The frozen expression-only PCA(50) plus k-means baseline was evaluated on the "
        "revised heterogeneous 50-gene/domain program at seven nominal fold changes "
        "for three independent SimSpace layout/count seeds. "
        "Within a seed, all random draws are paired across fold changes; only the "
        "domain-program multiplier changes.",
        "",
        f"First tested fold change with ARI >= 0.9: {thresholds}.",
        "",
        "Revised tier check: hard FC 2.0 median ARI "
        f"{tier_ari.loc[2.0]:.3f}; moderate FC 3.0 median ARI "
        f"{tier_ari.loc[3.0]:.3f} (seed range {moderate.min():.3f}--"
        f"{moderate.max():.3f}); easy FC 4.0 median ARI "
        f"{tier_ari.loc[4.0]:.3f}.",
        "",
        "## Adjusted Rand index",
        "",
        metric_table(frame, "adjusted_rand_index"),
        "",
        "## Normalized mutual information",
        "",
        metric_table(frame, "normalized_mutual_information"),
        "",
        "This is an engineering sweep over three seeds. The tested grid locates the "
        "transition but does not establish a universal detection threshold. The "
        "heterogeneous model spreads recovery across seeds and introduces a partial "
        "moderate regime, but hard PCA/k-means assignments can still change abruptly "
        "within an individual seed.",
        "",
    ]
    output.write_text("\n".join(lines))


def main() -> None:
    output_root = SCRIPT_DIR / "pilot_results"
    output_root.mkdir(parents=True, exist_ok=True)
    frame = run_sweep()
    if len(frame) != len(FOLD_CHANGES) * len(LAYOUT_SEEDS):
        raise RuntimeError("Incomplete expression signal sweep")
    if not frame["produced_cluster_count"].eq(4).all():
        raise RuntimeError("Expression baseline did not produce four clusters")
    frame.to_csv(
        output_root / "revised_expression_signal_sweep.tsv", sep="\t", index=False
    )
    render_sweep(frame, output_root / "revised_expression_signal_sweep.png")
    write_report(frame, output_root / "REVISED_EXPRESSION_SIGNAL_SWEEP.md")
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
