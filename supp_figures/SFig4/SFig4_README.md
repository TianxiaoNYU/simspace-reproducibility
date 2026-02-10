# Supplementary Figure 4

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 4**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig4.py`  
  Entry script to reproduce Supp Figure 4 (generates the figure and/or panel outputs).

- `Panel_C_D_data/`  
  Inputs and/or cached intermediates required for **Panel C, D**.
---

## Environment

### Default (SimSpace)
Most of Supp Figure 4 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Supp Figure 4

From the repository root:
```bash
python supp_figures/SFig4/SFig4.py
```

Or from this directory:
```bash
python SFig4.py
```

---

## Panel C D data

### BANKSY_xenium_domain.csv

They records the spatial clustering results for the Xenium reference and SimSpace data. Details can be found in `main_figures/Fig3/`

### Xenium_reference_count.csv
External data collected from High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis [https://www.nature.com/articles/s41467-023-43458-x]. This tile will used as the reference for scCube as it requires a reference to generate the molecular data even for its refernece-free mode.

### scCube_fitted_count.csv & Xenium_Breast_scCube.csv

They records the reference-based simulation results from scCube. Details can be found in `main_figures/Fig3/`

### simspace_fitted_count.csv
They records the reference-based simulation results from simspace reference-based simulations. Details can be found in `main_figures/Fig3/`

---
