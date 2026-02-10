# Supplementary Figure 7

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 7**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig7.py`  
  Entry script to reproduce Supp Figure 7 (generates the figure and/or panel outputs).

- `data/`  
  Inputs and/or cached intermediates required for panels.
---

## Environment

### Default (SimSpace)
Most of Supp Figure 7 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Supp Figure 7

From the repository root:
```bash
python supp_figures/SFig7/SFig7.py
```

Or from this directory:
```bash
python SFig7.py
```

---

## data

### fitted_pvals.csv

Cache file which stores the p-values or similar statistics from all 7 SVG detection methods. It will be used to generate the heatmap in panel B.

### METHOD_NAME/\*_280.csv, METHOD_NAME/\*_281.csv, ...

For all 7 SVG detection methods, we stored their SVG results for all 8 synthetic datasets we tested in the reference-free mode. These files contain the gene ID and the computed SVG statistics from each method.

### METHOD_NAME/\*_ref.csv, METHOD_NAME/\*_xenium.csv, ...

These files contain the SCG detection results in the reference-based mode, where `*_ref.csv` records the SVG results on the SimSpace simulated data, and `*_xenium.csv` records the SVG results on Xenium reference. 
---
