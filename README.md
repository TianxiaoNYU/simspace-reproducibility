# SimSpace Reproducibility Repository

This repository contains scripts and data to reproduce **all main and supplementary figures** for the SimSpace manuscript.

- Main figures: `main_figures/`
- Supplementary figures: `supp_figures/`
- Example outputs: `example_output/`
- Figure-to-script mapping: `FIGURE_INDEX.md`
- Default environment: `environment.yml`

---

## Quick start

### 1) Create the default environment
Using conda:
```bash
conda env create -f environment.yml
conda activate simspace-repro
pip install simspace==0.3.1
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
Running under: macOS Sequoia 15.7.3

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
 [1] splatter_1.30.0             SRTsim_0.99.8               SingleCellExperiment_1.28.1
 [4] SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0       
 [7] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
[10] BiocGenerics_0.52.0         MatrixGenerics_1.18.1       matrixStats_1.5.0          

loaded via a namespace (and not attached):
  [1] DBI_1.2.3               deldir_2.0-4            gridExtra_2.3           rlang_1.1.6            
  [5] magrittr_2.0.3          shinydashboard_0.7.3    spatstat.geom_3.4-1     e1071_1.7-16           
  [9] compiler_4.4.2          vctrs_0.6.5             pkgconfig_2.0.3         crayon_1.5.3           
 [13] fastmap_1.2.0           backports_1.5.0         magick_2.8.7            XVector_0.46.0         
 [17] promises_1.3.3          UCSC.utils_1.2.0        purrr_1.0.4             xfun_0.52              
 [21] zlibbioc_1.52.0         Rvcg_0.25               jsonlite_2.0.0          later_1.4.2            
 [25] DelayedArray_0.32.0     BiocParallel_1.40.2     spatstat.utils_3.1-4    pdist_1.2.1            
 [29] broom_1.0.8             parallel_4.4.2          R6_2.6.1                spatstat.data_3.1-6    
 [33] RColorBrewer_1.1-3      bezier_1.1.2            spatstat.univar_3.1-3   car_3.1-3              
 [37] Rcpp_1.0.14             iterators_1.0.14        knitr_1.50              base64enc_0.1-3        
 [41] FNN_1.1.4.1             httpuv_1.6.16           Matrix_1.7-4            tidyselect_1.2.1       
 [45] viridis_0.6.5           rstudioapi_0.17.1       dichromat_2.0-0.1       abind_1.4-8            
 [49] spatstat.random_3.4-1   doParallel_1.0.17       codetools_0.2-20        lattice_0.22-7         
 [53] tibble_3.3.0            shiny_1.10.0            S7_0.2.0                evaluate_1.0.3         
 [57] sf_1.0-21               polyclip_1.10-7         units_0.8-7             proxy_0.4-27           
 [61] BiocManager_1.30.26     pillar_1.10.2           ggpubr_0.6.0            carData_3.0-5          
 [65] KernSmooth_2.23-26      checkmate_2.3.2         DT_0.33                 foreach_1.5.2          
 [69] Morpho_2.12             plotly_4.10.4           generics_0.1.4          sp_2.2-0               
 [73] colorRamps_2.3.4        ggplot2_4.0.1           scales_1.4.0            xtable_1.8-4           
 [77] class_7.3-23            glue_1.8.0              lazyeval_0.2.2          tools_4.4.2            
 [81] data.table_1.17.4       locfit_1.5-9.12         ggsignif_0.6.4          rgl_1.3.18             
 [85] grid_4.4.2              tidyr_1.3.1             shinyBS_0.61.1          colorspace_2.1-1       
 [89] GenomeInfoDbData_1.2.13 Formula_1.2-5           cli_3.6.5               S4Arrays_1.6.0         
 [93] viridisLite_0.4.2       dplyr_1.1.4             concaveman_1.1.0        gtable_0.3.6           
 [97] rstatix_0.7.2           digest_0.6.37           classInt_0.4-11         SparseArray_1.6.2      
[101] htmlwidgets_1.6.4       farver_2.1.2            htmltools_0.5.8.1       lifecycle_1.0.4        
[105] httr_1.4.7              mime_0.13               MASS_7.3-65            
```