# Figure 4

This folder contains everything needed to reproduce the code-generated panels in **Figure 4**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `Fig4.py`  
  Entry script to reproduce Figure 4 (generates the figure and/or panel outputs).

- `Panel_A_data/`  
  Inputs and/or cached intermediates required for **Panel A**.

- `Panel_B_C_D_data/`
  Inputs and/or cached intermediates required for **Panel B, C, D**.  

- `Panel_B_C_D_src/`
  External code for Panel B, C, D that requires a separate environment for all deconvolution methods.  
---

## Environment

### Default (SimSpace)
Most of Figure 4 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

---

## Run Figure 4

From the repository root:
```bash
python main_figures/Fig4/Fig4.py
```

Or from this directory:
```bash
python Fig4.py
```

---

## Panel A data 

### simspace_fitted_params.json

SimSpace parameter file. Used for generating example plots in the workflow

---

## Panel B C D data 

### jaccard

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the jaccard index for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 


### pcc

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the pearson correlations for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 


### rmse

kernel{i}_niche{j}_state{k}_seed{l}.csv: 
- the result table containing the rooted mean square error for every spot's deconvolution results in the simulated data with {i, j, k, l} parameters. 

### Xenium_reference_*.csv

External data collected from High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis [https://www.nature.com/articles/s41467-023-43458-x]. This tile will used as the reference for scCube as it requires a reference to generate the molecular data even for its refernece-free mode.

Different files represents tiles with different number of cell types for SimSpace's simulation and the downstream deconvolution.

---

## Panel B C D src

### simspace_convolve.py

The script to generate SimSpace simulation data for the deconvolution benchmarking. One can run that using `simspace-repro` conda env. It will generate all the convolved spot-level spatial data (examples in `tmp_convolve_data`) for the downsteam deconvolution algorithms.

### tmp_convolve_data/meta and tmp_convolve_data/omics

SimSpace's simulated spot-level metadata and expression data. Stored in the format as `spot_meta_kernel{i}_niche{j}_state{k}_seed{l}.csv` and `spot_omics_kernel{i}_niche{j}_state{k}_seed{l}.csv`, where {i, j, k, l} represents the parameters used while simulating.

### deconvolution/main.py

The key script to (a) deconvolve the simulated convolved data using 6 differenet methods (b) compute the results statistics (Jaccard index, PCC, RMSE). To run the whole script, one should have all 6 methods, CARD, cell2location, RCTD, Seurat, spatialDWLS, and STdeconvolve installed in their own environment. 

### deconvolution/*.R (or *.py)

The individual scripts to generate the corresponding deconvolution results