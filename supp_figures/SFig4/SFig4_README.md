# Supplementary Figure 4

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 4**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig4.py`  
  Entry script to reproduce Supp Figure 4 (generates the figure and/or panel outputs).

- `Panel_A_B_C_data/`  
  Inputs and/or cached intermediates required for **Panel A, B, C**.

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

## Panel A B C data

### MERFISH
External data downloaded from the study: Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region [https://www.science.org/doi/10.1126/science.aau5324].

### Xenium
External data downloaded from the study: High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis [https://www.nature.com/articles/s41467-023-43458-x].
