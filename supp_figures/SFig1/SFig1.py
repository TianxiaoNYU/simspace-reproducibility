#!/usr/bin/env python3
"""
Simple scaling plots from benchmark results.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse, os

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default='results.csv')
    p.add_argument("--outdir", default="./")
    return p.parse_args()

def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    df = pd.read_csv(os.path.join(script_dir, args.csv))
    os.makedirs(args.outdir, exist_ok=True)
    

    agg = df.groupby(["method","grid_n"]).agg(
        runtime_mean=("runtime_s","mean"),
        runtime_std=("runtime_s","std"),
        mem_mean=("peak_ram_mb","mean"),
        mem_std=("peak_ram_mb","std")
    ).reset_index()

    # Runtime plot
    plt.figure(figsize=(5,4), dpi=250)
    for m in agg["method"].unique():
        d = agg[agg["method"]==m]
        plt.errorbar(d["grid_n"]**2, d["runtime_mean"], yerr=d["runtime_std"], marker='o', label=m)
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("# cells (grid size^2)"); plt.ylabel("Runtime (s)")
    plt.legend(); plt.tight_layout()
    plt.savefig(f"{args.outdir}/SFig1_panel_B1.png")

    # Memory plot
    plt.figure(figsize=(5,4), dpi=250)
    for m in agg["method"].unique():
        d = agg[agg["method"]==m]
        plt.errorbar(d["grid_n"]**2, d["mem_mean"], yerr=d["mem_std"], marker='o', label=m)
    plt.xscale("log")
    plt.xlabel("# cells (grid size^2)"); plt.ylabel("Peak RAM (MB)")
    plt.legend(); plt.tight_layout()
    plt.savefig(f"{args.outdir}/SFig1_panel_B2.png")

if __name__ == "__main__":
    main()