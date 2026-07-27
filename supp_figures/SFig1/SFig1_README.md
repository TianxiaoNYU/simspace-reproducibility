# Supplementary Figure 1

This folder contains everything needed to reproduce the code-generated panels in **Supp Figure 1**, including panel-specific data and (when needed) panel-specific external code.

---

## Contents

- `SFig1.py`  
  Entry script to reproduce Supp Figure 1 (generates the figure and/or panel outputs).

- `data/`
  Tracked 220-gene scCube source matrices and matching cell metadata.

- `generate_sccube_runtime_inputs.py`
  Deterministically expands the tracked scCube matrices to 2,200
  feature rows for the runtime stress test.

- `generated_data_manifest.tsv`
  Frozen generation parameters, dimensions, seeds, NumPy version, file
  sizes, and SHA-256 checksums for the generated scCube matrices.

- `generated_data/`
  Locally generated 2,200-feature scCube inputs. This large directory is
  intentionally excluded from Git.

- `runners.py`
  Python runners used to obtain runtime and memory benchmark results.

- `results*.csv`
  Simulation runtime and memory consumption for all three methods under
  different simulation sizes.
---

## Generate the 2,200-feature scCube runtime inputs

The expanded matrices are synthetic dimensionality stress-test inputs
for runtime benchmarking. They are not intended to represent a
biologically realistic 2,200-gene expression panel.

Each original profile is retained exactly. Nine nonnegative variants
are generated using deterministic Gaussian jitter with standard
deviation
`max(0.05 * source-gene SD, 0.001)`. File-specific seeds are `2200 +
grid_size`. The script also copies the matching cell metadata and writes
the checksum manifest.

From this directory:

```bash
python generate_sccube_runtime_inputs.py
```

The generator refuses to replace existing outputs by default. To verify
an existing `generated_data/` directory exactly without rewriting it:

```bash
python generate_sccube_runtime_inputs.py --verify-existing
```

To deliberately regenerate the directory:

```bash
python generate_sccube_runtime_inputs.py --overwrite
```

Before running the scCube benchmark, generate or verify these inputs.
Do not add `generated_data/` to Git; commit the generator, tracked
220-gene inputs, manifest, benchmark results, and documentation instead.

## Environment

### Default (SimSpace)
Most of Supp Figure 1 should run in the default environment:
```bash
# Run the line below if not creating the repro env before
# conda env create -f ../../environment.yml
conda activate simspace-repro
```

The environment pins `memory_profiler=0.61.0`, which is required for
peak-RAM measurement by the benchmark runner.

---

## Run Supp Figure 1

By default, the plotting script reads `results.csv` for the 220-gene
benchmark and `results_2200_gene.csv` for the 2,200-gene benchmark. It
writes the following matched panels to this directory:

- `SFig1_panel_B1.png`: runtime, 220 genes
- `SFig1_panel_B2.png`: peak RAM, 220 genes
- `SFig1_panel_C1.png`: runtime, 2,200 genes
- `SFig1_panel_C2.png`: peak RAM, 2,200 genes

From the repository root:
```bash
python supp_figures/SFig1/SFig1.py
```

Or from this directory:
```bash
python SFig1.py
```

---

### Benchmark simspace
```bash
conda activate simspace-repro
python bench.py \
  --methods simspace \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 220 \
  --outfile results_simspace.csv
```

### Benchmark sccube
```bash
# Run the line below if not creating the sccube env before
# conda env create -f sccube_environment.yml
conda activate sccube
python generate_sccube_runtime_inputs.py --verify-existing
python bench.py \
  --methods sccube \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 220 \
  --sccube-data-dir data \
  --outfile results_sccube.csv
```

### Benchmark scMultiSim

Install scMultiSim in the local R environment before running this benchmark.

```bash
conda activate simspace-repro
conda activate sccube
python bench.py \
  --methods scmultisim \
  --grid-sizes 10 15 20 25 \
  --repeats 5 \
  --threads 1 \
  --n-genes 220 \
  --outfile results_scmultisim.csv
```

### Benchmark the 2,200-gene condition

Use the same environments and grid sizes shown above for SimSpace and
scCube, replacing the gene count and output names:

```bash
python bench.py \
  --methods simspace \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --outfile results_simspace_2200.csv

python generate_sccube_runtime_inputs.py --verify-existing
python bench.py \
  --methods sccube \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --sccube-data-dir generated_data \
  --outfile results_sccube_2200.csv

python bench.py \
  --methods scmultisim \
  --grid-sizes 10 12 15 \
  --repeats 5 \
  --threads 1 \
  --n-genes 2200 \
  --outfile results_scmultisim_2200.csv
```

The planned 18 x 18 scMultiSim run was stopped manually because the
15 x 15 condition already required approximately 500 seconds per run;
no 18 x 18 result is included.

Below is the session info for scMultiSim:
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
[1] scMultiSim_1.2.0

loaded via a namespace (and not attached):
 [1] generics_0.1.4              SparseArray_1.6.2           lattice_0.22-7             
 [4] digest_0.6.37               magrittr_2.0.3              grid_4.4.2                 
 [7] RColorBrewer_1.1-3          iterators_1.0.14            foreach_1.5.2              
[10] jsonlite_2.0.0              Matrix_1.7-4                ape_5.8-1                  
[13] GenomeInfoDb_1.42.3         httr_1.4.7                  UCSC.utils_1.2.0           
[16] scales_1.4.0                codetools_0.2-20            abind_1.4-8                
[19] cli_3.6.5                   rlang_1.1.6                 crayon_1.5.3               
[22] XVector_0.46.0              Biobase_2.66.0              DelayedArray_0.32.0        
[25] Rtsne_0.17                  S4Arrays_1.6.0              tools_4.4.2                
[28] parallel_4.4.2              BiocParallel_1.40.2         dplyr_1.1.4                
[31] zeallot_0.2.0               ggplot2_4.0.1               GenomeInfoDbData_1.2.13    
[34] SummarizedExperiment_1.36.0 BiocGenerics_0.52.0         vctrs_0.6.5                
[37] R6_2.6.1                    matrixStats_1.5.0           stats4_4.4.2               
[40] lifecycle_1.0.4             zlibbioc_1.52.0             S4Vectors_0.44.0           
[43] IRanges_2.40.1              pkgconfig_2.0.3             pillar_1.10.2              
[46] gtable_0.3.6                glue_1.8.0                  Rcpp_1.0.14                
[49] tibble_3.3.0                GenomicRanges_1.58.0        tidyselect_1.2.1           
[52] rstudioapi_0.17.1           MatrixGenerics_1.18.1       dichromat_2.0-0.1          
[55] farver_2.1.2                nlme_3.1-168                compiler_4.4.2             
[58] S7_0.2.0                    markdown_2.0 
```
