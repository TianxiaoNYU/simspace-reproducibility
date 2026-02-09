# Environment Information

### BANKSY session info
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
 [1] cowplot_1.1.3               scater_1.34.1               ggplot2_4.0.1              
 [4] scuttle_1.16.0              SpatialExperiment_1.16.0    SingleCellExperiment_1.28.1
 [7] SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0       
[10] GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
[13] BiocGenerics_0.52.0         MatrixGenerics_1.18.1       matrixStats_1.5.0          
[16] Banksy_1.2.0               

loaded via a namespace (and not attached):
 [1] beeswarm_0.4.0          gtable_0.3.6            rjson_0.2.23            ggrepel_0.9.6          
 [5] lattice_0.22-7          vctrs_0.6.5             tools_4.4.2             generics_0.1.4         
 [9] parallel_4.4.2          tibble_3.3.0            sccore_1.0.6            pkgconfig_2.0.3        
[13] BiocNeighbors_2.0.1     Matrix_1.7-4            data.table_1.17.4       RColorBrewer_1.1-3     
[17] S7_0.2.0                lifecycle_1.0.4         GenomeInfoDbData_1.2.13 compiler_4.4.2         
[21] farver_2.1.2            aricode_1.0.3           codetools_0.2-20        vipor_0.4.7            
[25] pillar_1.10.2           crayon_1.5.3            BiocParallel_1.40.2     uwot_0.2.3             
[29] DelayedArray_0.32.0     dbscan_1.2.2            viridis_0.6.5           magick_2.8.7           
[33] abind_1.4-8             mclust_6.1.1            rsvd_1.0.5              tidyselect_1.2.1       
[37] BiocSingular_1.22.0     dplyr_1.1.4             grid_4.4.2              cli_3.6.5              
[41] SparseArray_1.6.2       magrittr_2.0.3          leidenAlg_1.1.5         S4Arrays_1.6.0         
[45] dichromat_2.0-0.1       withr_3.0.2             scales_1.4.0            UCSC.utils_1.2.0       
[49] ggbeeswarm_0.7.2        XVector_0.46.0          httr_1.4.7              igraph_2.1.4           
[53] gridExtra_2.3           ScaledMatrix_1.14.0     beachmat_2.22.0         RcppHungarian_0.3      
[57] viridisLite_0.4.2       irlba_2.3.5.1           rlang_1.1.6             Rcpp_1.0.14            
[61] glue_1.8.0              rstudioapi_0.17.1       jsonlite_2.0.0          R6_2.6.1               
[65] zlibbioc_1.52.0 
```

### SEDR conda env info
- We followed the environment requirements from the GitHub official page of SEDR, which has the following pakcages:
```bash
python==3.11.3
torch==2.0.1
cudnn==12.1
numpy==1.24.3
scanpy==1.9.3
anndata==0.9.1
ryp2==3.5.12
pandas==2.0.1
scipy==1.10.1
scikit-learn==1.2.2
tqdm==4.65.0
matplotlib==3.7.1
seaborn==0.12.2
R==4.0.3
```

### SpatialPCA session info
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4      tidyr_1.3.1      ggplot2_4.0.1    SpatialPCA_1.3.0

loaded via a namespace (and not attached):
  [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6           
  [5] magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.4-1    matrixStats_1.5.0     
  [9] ggridges_0.5.6         compiler_4.4.2         png_0.1-8              vctrs_0.6.5           
 [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0         
 [17] promises_1.3.3         pracma_2.4.4           purrr_1.0.4            jsonlite_2.0.0        
 [21] goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4   pdist_1.2.1           
 [25] irlba_2.3.5.1          parallel_4.4.2         cluster_2.1.8.1        R6_2.6.1              
 [29] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-6   
 [33] reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40         
 [37] scattermore_1.2        Rcpp_1.0.14            iterators_1.0.14       tensor_1.5            
 [41] future.apply_1.20.0    zoo_1.8-14             sctransform_0.4.2      matlab_1.0.4.1        
 [45] httpuv_1.6.16          Matrix_1.7-4           splines_4.4.2          igraph_2.1.4          
 [49] tidyselect_1.2.1       abind_1.4-8            rstudioapi_0.17.1      dichromat_2.0-0.1     
 [53] spatstat.random_3.4-1  doParallel_1.0.17      codetools_0.2-20       miniUI_0.1.2          
 [57] spatstat.explore_3.4-3 listenv_0.9.1          lattice_0.22-7         tibble_3.3.0          
 [61] plyr_1.8.9             withr_3.0.2            shiny_1.10.0           S7_0.2.0              
 [65] askpass_1.2.1          ROCR_1.0-11            Rtsne_0.17             future_1.58.0         
 [69] fastDummies_1.7.5      survival_3.8-3         CompQuadForm_1.4.3     polyclip_1.10-7       
 [73] fitdistrplus_1.2-2     pillar_1.10.2          Seurat_5.4.0           KernSmooth_2.23-26    
 [77] foreach_1.5.2          plotly_4.10.4          generics_0.1.4         RcppHNSW_0.6.0        
 [81] sp_2.2-0               scales_1.4.0           globals_0.18.0         xtable_1.8-4          
 [85] glue_1.8.0             lazyeval_0.2.2         tools_4.4.2            data.table_1.17.4     
 [89] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.1.3         
 [93] grid_4.4.2             umap_0.2.10.0          nlme_3.1-168           patchwork_1.3.0       
 [97] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            SPARK_1.1.1           
[101] viridisLite_0.4.2      uwot_0.2.3             gtable_0.3.6           digest_0.6.37         
[105] progressr_0.15.1       ggrepel_0.9.6          htmlwidgets_1.6.4      SeuratObject_5.1.0    
[109] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
[113] mime_0.13              openssl_2.3.3          MASS_7.3-65    
```