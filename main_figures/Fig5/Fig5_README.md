# Figure 5

This folder contains everything needed to reproduce the code-generated panels in **Figure 5**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `Fig5.py`  
  Entry script to reproduce Figure 5 (generates the figure and/or panel outputs).

- `Panel_B_C_D_data/`  
  Inputs and/or cached intermediates required for **Panel B, C, D**.

- `Panel_B_C_D_src/`
  External code for Panel B, C, D that requires a separate environment for the SVG detection methods.  
---

## Environment

### Default (SimSpace)
Most of Figure 5 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Figure 5

From the repository root:
```bash
python main_figures/Fig5/Fig5.py
```

Or from this directory:
```bash
python Fig5.py
```

---

## Panel B C D data 

### fitted_pvals.csv

Cache file which stores the p-values or similar statistics from all 7 SVG detection methods. It will be used to generate the heatmap in panel B.

### METHOD_NAME/\*_280.csv, METHOD_NAME/\*_281.csv, ...

For all 7 SVG detection methods, we stored their SVG results for all 8 synthetic datasets we tested in the reference-free mode. These files contain the gene ID and the computed SVG statistics from each method.

### METHOD_NAME/\*_ref.csv, METHOD_NAME/\*_xenium.csv, ...

These files contain the SCG detection results in the reference-based mode, where `*_ref.csv` records the SVG results on the SimSpace simulated data, and `*_xenium.csv` records the SVG results on Xenium reference. 

---

## Panel B C D src

### Xenium_reference_count.csv & Xenium_reference_metadata.csv
External data collected from High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis [https://www.nature.com/articles/s41467-023-43458-x]. This tile will used as the reference for SimSpace to do the reference-based simulations.

### simspace_data_prep.py

Generate the reference-free simulation datasets, which will be stored in `benchmark_datasets` for the downstream analysis.

### benchmark_datasets

Store all the SimSpace generated data for SVG detection comparison. SimSpace's reference-based simulation results are within `./Xenium_fitted/`, and the Xenium reference within `./Xenium_red/`

### panel_B.R

R script to draw the heatmap in panel B of Figure 5, using the cache `Panel_B_C_D_data/fitted_pvals.csv`. R packages `ggplot2` and `pheatmap` are necessary.

### METHOD_NAME/

The script of each deconvolution methods, including scBSP, Giotto, Hotspot, Moran's I, SPARK, and spatialDE. We followed the provided guidance on each package's official page to have them installed, so it is recommended to follow the same guidance to reproduce the results if one desired.