#!/usr/bin/env python3
"""Generate the 220- and 2,200-gene runtime benchmark panels."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


METHOD_STYLES = {
    "sccube": {
        "label": "scCube",
        "color": "#0072B2",
        "marker": "o",
    },
    "scmultisim": {
        "label": "scMultiSim",
        "color": "#E69F00",
        "marker": "s",
    },
    "simspace": {
        "label": "SimSpace",
        "color": "#009E73",
        "marker": "^",
    },
}
METHOD_ORDER = tuple(METHOD_STYLES)
REQUIRED_COLUMNS = {
    "method",
    "grid_n",
    "runtime_s",
    "peak_ram_mb",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate matched runtime and peak-memory panels for the "
            "220- and 2,200-gene benchmarks."
        )
    )
    p.add_argument(
        "--csv-220",
        default="results.csv",
        help="CSV containing the 220-gene benchmark results.",
    )
    p.add_argument(
        "--csv-2200",
        default="results_2200_gene.csv",
        help="CSV containing the 2,200-gene benchmark results.",
    )
    p.add_argument(
        "--outdir",
        default=None,
        help=(
            "Directory in which to write panels B1, B2, C1, and C2 "
            "(default: this script's directory)."
        ),
    )
    return p.parse_args()


def resolve_input_path(script_dir, value):
    path = Path(value)
    return path if path.is_absolute() else script_dir / path


def load_and_aggregate(path):
    df = pd.read_csv(path)
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"{path} is missing required columns: {missing_list}")

    agg = (
        df.groupby(["method", "grid_n"])
        .agg(
            runtime_mean=("runtime_s", "mean"),
            runtime_std=("runtime_s", "std"),
            mem_mean=("peak_ram_mb", "mean"),
            mem_std=("peak_ram_mb", "std"),
        )
        .reset_index()
    )
    agg[["runtime_std", "mem_std"]] = agg[
        ["runtime_std", "mem_std"]
    ].fillna(0)
    return agg


def ordered_methods(*aggregates):
    observed = {
        method
        for aggregate in aggregates
        for method in aggregate["method"].unique()
    }
    known = [method for method in METHOD_ORDER if method in observed]
    return known + sorted(observed.difference(METHOD_ORDER))


def shared_limits(aggregates):
    all_results = pd.concat(aggregates, ignore_index=True)
    cell_counts = all_results["grid_n"] ** 2

    runtime_lower = (
        all_results["runtime_mean"] - all_results["runtime_std"]
    ).clip(lower=all_results["runtime_mean"] * 0.5)
    runtime_upper = (
        all_results["runtime_mean"] + all_results["runtime_std"]
    )
    memory_upper = all_results["mem_mean"] + all_results["mem_std"]

    return {
        "x": (cell_counts.min() / 1.2, cell_counts.max() * 1.2),
        "runtime": (runtime_lower.min() / 1.2, runtime_upper.max() * 1.2),
        "memory": (0, memory_upper.max() * 1.08),
    }


def plot_panel(
    aggregate,
    methods,
    metric,
    ylabel,
    gene_label,
    output_path,
    x_limits,
    y_limits,
    log_y=False,
):
    fig, ax = plt.subplots(figsize=(5, 4), dpi=250)

    for index, method in enumerate(methods):
        data = aggregate[aggregate["method"] == method].sort_values("grid_n")
        if data.empty:
            continue

        style = METHOD_STYLES.get(
            method,
            {
                "label": method,
                "color": f"C{index}",
                "marker": "o",
            },
        )
        ax.errorbar(
            data["grid_n"] ** 2,
            data[f"{metric}_mean"],
            yerr=data[f"{metric}_std"],
            color=style["color"],
            marker=style["marker"],
            markersize=6,
            linewidth=1.8,
            capsize=3,
            label=style["label"],
        )

    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_xlabel(r"Number of cells (grid size$^2$)")
    ax.set_ylabel(ylabel)
    ax.set_title(gene_label)
    ax.legend(frameon=False)
    ax.tick_params(direction="out")
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    outdir = Path(args.outdir) if args.outdir else script_dir
    outdir.mkdir(parents=True, exist_ok=True)

    results_220 = load_and_aggregate(
        resolve_input_path(script_dir, args.csv_220)
    )
    results_2200 = load_and_aggregate(
        resolve_input_path(script_dir, args.csv_2200)
    )
    methods = ordered_methods(results_220, results_2200)
    limits = shared_limits([results_220, results_2200])

    panel_specs = (
        (
            results_220,
            "runtime",
            "Runtime (s)",
            "220 genes",
            "SFig1_panel_B1.png",
            limits["runtime"],
            True,
        ),
        (
            results_220,
            "mem",
            "Peak RAM (MB)",
            "220 genes",
            "SFig1_panel_B2.png",
            limits["memory"],
            False,
        ),
        (
            results_2200,
            "runtime",
            "Runtime (s)",
            "2,200 genes",
            "SFig1_panel_C1.png",
            limits["runtime"],
            True,
        ),
        (
            results_2200,
            "mem",
            "Peak RAM (MB)",
            "2,200 genes",
            "SFig1_panel_C2.png",
            limits["memory"],
            False,
        ),
    )

    for (
        aggregate,
        metric,
        ylabel,
        gene_label,
        filename,
        y_limits,
        log_y,
    ) in panel_specs:
        plot_panel(
            aggregate=aggregate,
            methods=methods,
            metric=metric,
            ylabel=ylabel,
            gene_label=gene_label,
            output_path=outdir / filename,
            x_limits=limits["x"],
            y_limits=y_limits,
            log_y=log_y,
        )

    print(f"Wrote four benchmark panels to {outdir.resolve()}")


if __name__ == "__main__":
    main()
