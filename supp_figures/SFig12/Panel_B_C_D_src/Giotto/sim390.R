library(tidyverse)
library(SingleCellExperiment) #provides infrastructure for single-cell analysis
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
#library(zellkonverter)
library(readr)
library(dplyr)
library(tidyr)
library(Matrix)
library(Giotto)

library(here)
library(ggplot2)
library(pheatmap)
library(data.table)

library(patchwork)
library(viridisLite)
library(viridis)

library(scuttle)
library(scran)

source("call_variable_genes.R")


meta390 <- read_csv("dataset_0515/sim_meta_niche3_state9_seed0.csv")

omics390 <- read_csv("dataset_0515/sim_omics_niche3_state9_seed0.csv")

countMat390 <- as.matrix(omics390)

dim(countMat390)

rownames(countMat390) <- rownames(meta390)


spatial_locs <- data.frame(
  x = meta390$col,
  y = meta390$row,
  row.names = rownames(countMat390)
)

meta390 <- as.data.frame(meta390)

rownames(spatial_locs) <- rownames(countMat390)

countMat390 <- t(countMat390)


giotto_obj <- createGobject(
  countMat = countMat390,
  spatial_locs = spatial_locs,
  cell_metadata = meta390,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "sim/svg_giotto_kmeans390.csv")

fwrite(svg_giotto_rank, file = "sim/svg_giotto_rank390.csv")