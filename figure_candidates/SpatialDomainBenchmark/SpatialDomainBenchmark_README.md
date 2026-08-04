# Spatial-domain benchmark figure candidate

Prepared and tested on 2026-08-03 on macOS arm64. These environments are
separate from `simspace-repro` and do not install or modify SimSpace.

This directory is a temporary, self-contained figure workspace. The benchmark
has not yet been assigned a main- or supplementary-figure number; when that is
decided, the directory can move unchanged to `main_figures/FigX` or
`supp_figures/SFigX` and adopt the corresponding figure filenames.

## Directory contents

- `environment.yml`: GraphST, STAGATE, SpaGCN, and spCLUE
- `environment-stlearn.yml`: isolated current stLearn runtime
- `install-methods.sh`: immutable source installs plus the host-R bridge
- `smoke-test.py` and `smoke-test-stlearn.py`: minimal runtime validation
- `resolved-osx-arm64*.yml`: exact lock-style snapshots from validation
- `BENCHMARK_DESIGN.md`: scientific, fairness, execution, and figure design
- `benchmark-config.draft.yml`: machine-readable production design
- `generate_pilot.py`: shared reference-free simulation implementation with
  four domains and ten cell types
- `run_*.py` and `run_banksy.R`: truth-blind method/control adapters and orchestrator
- `evaluate_pilot.py`: held-out ARI, NMI, boundary, runtime, and memory summary
- `pilot_data/` and `pilot_results/`: generated local engineering-pilot artifacts
- `generate_experiment.py`: production generator for the 18 paired datasets
- `run_experiment.py`: production orchestrator for all 108 named-method runs
- `evaluate_experiment.py`: held-out production aggregation and validation
- `experiment_data/` and `experiment_results/`: generated production artifacts

Each method receives only counts, generic gene identifiers, and coordinates.
Truth labels, the full gene-class table, and the latent molecular mean remain
in held-out provenance files read only by the evaluation/generation checks.

## Why there are two environments

Two is the minimum compatible environment count. The current stable
[stLearn 1.4.1 source release](https://github.com/BiomedicalMachineLearning/stLearn/tree/V1.4.1)
requires Python >=3.12, NumPy >=2.4, Scanpy >=1.12, SciPy >=1.17, and
scikit-learn >=1.8. The other four methods use the Scanpy 1.9 / AnnData 0.8
generation; for example, [GraphST](https://github.com/JinmiaoChenLab/GraphST)
documents Python 3.8 and NumPy 1.22.3, while
[spCLUE](https://github.com/EnchantedJoy/spCLUE) documents Python 3.9 and
NumPy 1.23.5. Those requirements cannot be solved faithfully in one runtime.

The integrated four-method environment uses spCLUE's scientific stack as the
baseline and Python 3.10 for current Apple-Silicon binary availability.
PyTorch 2.4.0 and PyG 2.6.1 are intentional compatibility lifts: this is the
first complete osx-arm64 combination in the PyG channel with working
`torch_sparse`, `torch_scatter`, and `torch_cluster` binaries. PyTorch is
channel-qualified because the compiled PyG extensions require the official
PyTorch-channel ABI.

The shared environment also pins `setuptools<81`: the pinned legacy
`umap-learn` imports `pkg_resources`, which setuptools 81 and later removed.

## Create and test

Create the four-method environment and then run its required source/host-R
post-install step:

```bash
mamba env create -f environment.yml
conda activate spatial-domain-benchmark
bash install-methods.sh
```

Create and test current stLearn separately:

```bash
mamba env create -f environment-stlearn.yml
conda activate spatial-domain-benchmark-stlearn
python smoke-test-stlearn.py
```

The two `resolved-osx-arm64*.yml` files are lock-style snapshots of the exact
versions installed during validation. Prefer the shorter primary YAMLs for a
new platform. The resolved four-method snapshot still needs
`install-methods.sh`, because rpy2 must bind to host R and upstream spCLUE
has no `setup.py` or `pyproject.toml`.

## Pinned methods

| Method | Version/source | Immutable revision |
|---|---|---|
| GraphST | 1.1.1 | `d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74` |
| STAGATE_pyG | 1.0.0 | `ae1158ca8cf1eb6bb8ee198298552d44c9ac21db` |
| SpaGCN | 1.2.7 | `dc7a1c26ea0fdf4dfe7064adc7699be141b4871f` |
| spCLUE | source only | `bbd2c342e7a67c1617275f721cec2e3f4c23a799` |
| stLearn | 1.4.1 release | `03e38844a599ef4f3df8235e3f57bd74a9e96fcf` |

Official installation sources: [GraphST](https://deepst-tutorials.readthedocs.io/en/latest/Installation.html),
[STAGATE_pyG](https://stagate.readthedocs.io/en/latest/Installation_pyG.html),
[SpaGCN](https://github.com/jianhuupenn/SpaGCN),
[spCLUE](https://github.com/EnchantedJoy/spCLUE), and
[stLearn](https://stlearn.readthedocs.io/en/stable/installation.html).

## R policy

No Conda R runtime is used. `install-methods.sh` checks the normal local R
library, installs `mclust` there only if missing, and builds rpy2 3.5.17 in ABI
mode against R found on `PATH`. On the validation host, the existing packages
were:

- R 4.4.2
- mclust 6.1.1
- BANKSY 1.2.0

BANKSY remains in the existing local R/SimSpace setup and is not duplicated in
either Python environment.

## Run the engineering pilot

From this directory, generate the MRF dataset, run all six named methods and the
two prespecified baselines, and aggregate the held-out metrics with:

```bash
python run_pilot.py
```

An interrupted run can reuse the frozen data and select methods explicitly:

```bash
python run_pilot.py --skip-generation --methods graphst banksy stlearn
```

A paired signal level can be regenerated or rerun explicitly, for example:

```bash
python run_pilot.py --difficulty-label moderate --domain-fold-change 2.5
```

The main outputs are `pilot_results/pilot_summary.tsv`,
`pilot_results/pilot_overview.png`, and `pilot_results/PILOT_REPORT.md`.
After all three paired signal pilots are available, `compare_signal_pilots.py`
also writes `pilot_results/pilot_signal_comparison.tsv`,
`pilot_results/pilot_signal_comparison.png`, and
`pilot_results/PILOT_SIGNAL_COMPARISON.md`.
The focused three-seed expression-baseline sweep can be reproduced with:

```bash
conda run --no-capture-output -n simspace-repro python sweep_expression_signal.py
```

The archived original-design sweep remains in
`pilot_results/expression_signal_sweep.*`. The revised heterogeneous-program
sweep writes `pilot_results/revised_expression_signal_sweep.tsv`,
`pilot_results/revised_expression_signal_sweep.png`, and
`pilot_results/REVISED_EXPRESSION_SIGNAL_SWEEP.md`.
This pilot validates the workflow and diagnoses simulation difficulty; it is
not the production benchmark and its truth metrics must not be used to tune
method parameters.

The MRF pilot uses a 68 x 68 nominal grid, retains exactly 3,000 of 4,624
locations (64.9%), and applies coordinate jitter with standard deviation 0.55
lattice units, matching the sampling treatment used for SFig12 panel A. Its
continuous irregular domains use a niche MRF with `phi=5.25`, 40 sweeps, and a
Manhattan radius-2 neighborhood. Cell-type composition assigns 40% probability
mass to each domain's enriched types and shares the remaining 60% across all
ten types. The revised molecular design uses 50 heterogeneous domain-program
genes per domain with domain-, gene-, and location-level log-effect weights.
Its production hard, moderate, and easy nominal fold changes are 2.0, 3.0, and
4.0.

## Production experiment (scripted, not run)

The frozen matrix contains two pattern types (`curved_layers` and
`irregular_mrf`), three difficulty levels (hard, moderate, and easy), three
layout seeds (0, 1, and 2), and six named methods. This produces 18 datasets,
108 named-method runs, and 18 expression-only negative-control runs. Both
patterns use the same revised molecular model,
including the 40% enriched/60% shared cell-type composition and heterogeneous
50-gene/domain expression program. Within each pattern and layout seed, all
latent draws are paired across difficulty and only the nominal domain fold
change varies.

Inspect every planned command without generating data or launching a method:

```bash
python run_experiment.py --dry-run
```

Launch the complete experiment from this directory only after the design is
approved:

```bash
python run_experiment.py
```

After an interruption, reuse the generated datasets and preserve every valid
completed result with:

```bash
python run_experiment.py --skip-generation --resume-successful
```

The orchestrator writes a frozen plan, per-run command and standard-stream
logs, explicit timeout/failure records, `experiment_summary.tsv`,
`evaluation_manifest.json`, and `EXPERIMENT_REPORT.md`. Production evaluation
accepts only 3,000 aligned assignments with exactly four nonempty clusters.
The expression-only PCA(50) plus k-means control is labeled separately from the
six named spatial methods and never contributes to their completion denominator.

## Validation performed

`smoke-test.py` imports all four shared-environment methods, executes a compiled
PyTorch-Geometric/`torch_sparse` operation, and loads local-R `mclust` through
rpy2. `smoke-test-stlearn.py` constructs a small stLearn object and completes a
two-domain k-means clustering. Both tests passed on the validation host.
