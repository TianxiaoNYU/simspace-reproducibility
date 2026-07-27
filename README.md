# SimSpace Reproducibility Repository

This repository contains scripts and data to reproduce **all main and supplementary figures** for the SimSpace manuscript.

- Python script for main figures: `plot_main.py`
- Python script for Supplementary figures: `plot_supp.py`
- Main figures: `main_figures/`
- Supplementary figures: `supp_figures/`
- Example outputs: `example_output/`
- Figure-to-script mapping: `FIGURE_INDEX.md`
- Default environment: `environment.yml`

---

## Quick start

### 1) Create the default environment
We recommend using conda to setup the environment:

```bash
git clone https://github.com/TianxiaoNYU/simspace-reproducibility.git
```

```bash
conda env create -f environment.yml
conda activate simspace-repro
pip install simspace==0.3.2
```

### 2) Setting Up the R Environment for Omics Simulation

SimSpace supports external omics profile simulation via R-based tools, including **scDesign3**, **SRTsim**, and **splatter**. These tools are optional but recommended if you want to simulate gene expression profiles in addition to spatial patterns.

To enable this functionality, please install the required R packages manually in your system R environment. A detailed R session info is attached at the end of the README.

Steps:
1.	Ensure that R (version 4.4 or compatible) is installed on your system. You can download it from CRAN.
2.	Open an R session and install the required packages:
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
devtools::install_github("xzhoulab/SRTsim")
```
```R
if (!require("devtools", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("splatter"))
```

### 3) Reproduce all main figures
```bash
python plot_main.py
```

### 4) Reproduce all supplementary figures
```bash
python plot_supp.py
```

---

## Reproducing a single figure

Each figure has its own folder and a single entry script.

Example (Figure 2):
```bash
python main_figures/Fig2/Fig2.py
```

See `FIGURE_INDEX.md` for the complete mapping.

---

## Repository layout

```
simspace-reproducibility/
  example_output/
    Fig1/
      Fig1_panel_A.png
      Fig1_panel_B1.png
      ...
    Fig2/
    ...
    SFig1/
      SFig1_panel_B1.png
      SFig1_panel_B2.png
      SFig1_panel_C1.png
      SFig1_panel_C2.png
    SFig2/
    ...
  main_figures/
    Fig1/
    Fig2/
      Fig2.py
      Fig2_README.md
      Panel_D_data/
      Panel_D_src/
    Fig3/
      Fig3.py
      Fig3_README.md
      Panel_*_data/
      Panel_*_src/
    ...
  supp_figures/
    SFig1/
    SFig2/
    ...
  plot_main.py
  plot_supp.py
  FIGURE_INDEX.md
  environment.yml
  README.md
  LICENSE
```

---

## Panel data vs. external code (`Panel_*_data` and `Panel_*_src`)

Some figures are organized by panel to keep dependencies and data clear.

- `Panel_*_data/`  
  Inputs and/or cached intermediates required to reproduce that panel.  
  These are intended to work with the **default** environment unless explicitly stated otherwise.

- `Panel_*_src/`  
  External code that is difficult to integrate into the default environment (e.g., torch/CUDA, specialized R packages).  
  Each `Panel_*_src/` directory should be self-contained and include:
  - a short `Fig*_README.md` with exact setup, run commands, and conda env/R session info
  - a single entrypoint script for each method/package where possible

Panels that require a separate environment are explicitly labeled in `FIGURE_INDEX.md` and in the figure-level `Fig*_README.md`.

---

## Outputs

By convention, scripts write outputs under the its own directory, and print the output paths upon completion. 

We also provided the expected results, storing in `example_output/`. One can compare the outputs with the example_output to ensure the reproducibility.

---

## Notes on determinism / reproducibility

We aim for reproducible figure regeneration under the provided environments.  
Some third-party tools may still exhibit minor run-to-run variation even with fixed seeds (e.g., parallelism, non-deterministic low-level libraries, or upstream package behavior).

To minimize variability:
- use the provided environment specs (`environment.yml` or per-panel env files)
- for methods/packages beyond SimSpace, follow the README in each figure folder
- keep seeds fixed where applicable
- prefer CPU execution if GPU nondeterminism is observed

Where relevant, we provide cached intermediate results / caches to avoid refitting and stabilize outputs.

---

## Troubleshooting

### Common issues
- **Package/version mismatch:** recreate the environment from the relevant env file.
- **Missing optional dependencies:** check the error message and install into the correct env.
- **Large files:** some intermediates may be excluded from git and generated/downloaded on first run.

### Getting help
If you encounter an issue reproducing a figure, please open a GitHub issue and include:
- figure name (e.g., Fig3, SFig2)
- the command you ran
- full traceback/error message
- OS + Python version


---

## Additional: scDesgin3, SRTsim, and splatter R session info

```
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS 26.5.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] splatter_1.30.0             SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
 [4] Biobase_2.66.0              GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
 [7] IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0        
[10] MatrixGenerics_1.18.1       matrixStats_1.5.0           SRTsim_0.99.8              
[13] scDesign3_1.5.0            

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_2.0.0          magrittr_2.0.3         
  [5] spatstat.utils_3.2-1    magick_2.8.7            farver_2.1.2            zlibbioc_1.52.0        
  [9] vctrs_0.6.5             base64enc_0.1-3         rstatix_0.7.2           htmltools_0.5.8.1      
 [13] S4Arrays_1.6.0          Morpho_2.12             broom_1.0.8             SparseArray_1.6.2      
 [17] Formula_1.2-5           KernSmooth_2.23-26      htmlwidgets_1.6.4       plotly_4.10.4          
 [21] mime_0.13               lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
 [25] Matrix_1.7-4            R6_2.6.1                fastmap_1.2.0           GenomeInfoDbData_1.2.13
 [29] shiny_1.10.0            digest_0.6.37           colorspace_2.1-1        bezier_1.1.2           
 [33] ggpubr_0.6.0            pdist_1.2.1             colorRamps_2.3.4        httr_1.4.7             
 [37] polyclip_1.10-7         abind_1.4-8             compiler_4.4.2          proxy_0.4-27           
 [41] doParallel_1.0.17       S7_0.2.0                backports_1.5.0         BiocParallel_1.40.2    
 [45] carData_3.0-5           viridis_0.6.5           DBI_1.2.3               ggsignif_0.6.4         
 [49] MASS_7.3-65             concaveman_1.1.0        DelayedArray_0.32.0     classInt_0.4-11        
 [53] tools_4.4.2             units_0.8-7             httpuv_1.6.16           glue_1.8.0             
 [57] nlme_3.1-168            promises_1.3.3          grid_4.4.2              sf_1.0-21              
 [61] checkmate_2.3.2         generics_0.1.4          gtable_0.3.6            spatstat.data_3.1-9    
 [65] Rvcg_0.25               shinyBS_0.61.1          class_7.3-23            tidyr_1.3.1            
 [69] data.table_1.17.4       sp_2.2-0                car_3.1-3               XVector_0.46.0         
 [73] spatstat.geom_3.7-0     foreach_1.5.2           pillar_1.10.2           later_1.4.2            
 [77] splines_4.4.2           dplyr_1.1.4             lattice_0.22-7          survival_3.8-3         
 [81] FNN_1.1.4.1             deldir_2.0-4            gamlss.data_6.0-6       tidyselect_1.2.1       
 [85] locfit_1.5-9.12         knitr_1.50              gridExtra_2.3           xfun_0.52              
 [89] shinydashboard_0.7.3    DT_0.33                 UCSC.utils_1.2.0        lazyeval_0.2.2         
 [93] evaluate_1.0.3          codetools_0.2-20        tibble_3.3.0            cli_3.6.5              
 [97] xtable_1.8-4            dichromat_2.0-0.1       Rcpp_1.0.14             spatstat.random_3.4-4  
[101] spatstat.univar_3.1-6   parallel_4.4.2          rgl_1.3.18              ggplot2_4.0.1          
[105] mclust_6.1.1            gamlss.dist_6.1-1       viridisLite_0.4.2       scales_1.4.0           
[109] e1071_1.7-16            gamlss_5.4-22           purrr_1.0.4             crayon_1.5.3           
[113] rlang_1.1.6
```
