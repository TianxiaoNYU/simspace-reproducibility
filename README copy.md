# SimSpace Reproducibility (Nature Communications submission)

This repository contains the **analysis and reproduction script** for the manuscript:

> **SimSpace: a comprehensive in-silico spatial omics data simulation and modeling framework**

It is intentionally separated from the core SimSpace package repository to keep the package lightweight while providing **full computational reproducibility** for all results and figures presented in the manuscript. Still, one can find the reproduction notebook in the main repo of simspace package, as well.

---

## What this repo provides

- ✅ Scripts / notebooks to reproduce **all main and supplementary figures**
- ✅ **Pinned software environments**, including **exact versions of competitor methods** used in benchmarking
- ✅ Clear **data provenance** for all inputs and intermediate files

**Core SimSpace package repo: https://github.com/TianxiaoNYU/simspace**

---

## Repository structure (high-level)

- `main_figures.ipynb` – reproduces all main manuscript figures
- `supp_figures.ipynb` – reproduces all Supplementary figures
- `scripts/reproduce_all.sh` – non-interactive run (generates all figures + logs)
- `FIGURE_INDEX.md` – mapping from **Figure/Panel → code location → output file**
- `data_provenance.md` – source and generation procedure for each input/intermediate
- `env/` – pinned environments
  - `environment.yml` (conda)
  - `renv.lock` (R) and `scripts/install_R_deps.R`
  - version snapshots: `python_freeze.txt`, `R_sessionInfo.txt`

Outputs are generated under `outputs/` (not tracked by git).

---

## Quick start (recommended)

### 1) Create environment for python (Conda)
```bash
conda env create -f env/environment.yml
conda activate simspace-repro
```

### 2) Download packages in R

Since some of the external tools used in the paper are wrapped in R, it is recommended to manually download these packages in user's own local environment. These packages and their versions are:

- Spatial simulation
    - scDesign3: 1.5.0
    - scMultiSim: 1.2.0
    - SRTsim: 0.99.8

- Spatial clustering
    - BANKSY
    - SEDR
    - SpatialPCA

- Cell type deconvolution
    - CARD
    - RCTD
    - Seurat
    - spatialDWLS
    - STdeconvolve

- SVG detection
    - scBSP
    - Hotspot
    - SPARK

### 3) Reproduce figures
After installing the environment and downloading data:
```bash
jupyter lab
```
Then open:
- ./main_figures.ipynb
- ./supp_figures.ipynb

Each notebook is designed to run top-to-bottom without manual edits.