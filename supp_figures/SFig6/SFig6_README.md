# Supplementary Figure 6

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 6**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig6.py`  
  Entry script to reproduce Supp Figure 6 (generates the figure and/or panel outputs).

- `Panel_A_B_C_data/`  
  Inputs and/or cached intermediates required for **Panel A, B, C**.
---

## Environment

### Default (SimSpace)
Most of Supp Figure 6 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Supp Figure 6

From the repository root:
```bash
python supp_figures/SFig6/SFig6.py
```

Or from this directory:
```bash
python SFig6.py
```

---

## Panel A B C data

### jaccard

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the jaccard index for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 


### pcc

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the pearson correlations for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 


### rmse

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the rooted mean square error for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 

---
