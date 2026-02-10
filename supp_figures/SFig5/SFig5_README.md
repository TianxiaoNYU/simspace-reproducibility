# Supplementary Figure 5

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 5**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig5.py`  
  Entry script to reproduce Supp Figure 5 (generates the figure and/or panel outputs).

- `Panel_A_B_C_data/`  
  Inputs and/or cached intermediates required for **Panel A, B, C**.
---

## Environment

### Default (SimSpace)
Most of Supp Figure 5 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Supp Figure 5

From the repository root:
```bash
python supp_figures/SFig5/SFig5.py
```

Or from this directory:
```bash
python SFig5.py
```

---

## Panel A B C data

### Xenium_reference_*.csv

External data collected from High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis [https://www.nature.com/articles/s41467-023-43458-x]. This tile will used as the reference for scCube as it requires a reference to generate the molecular data even for its refernece-free mode.

Different files represents tiles with different number of cell types for SimSpace's simulation and the downstream deconvolution.


---
