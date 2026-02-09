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


meta380 <- read_csv("dataset_0515/sim_meta_niche3_state8_seed0.csv")

omics380 <- read_csv("dataset_0515/sim_omics_niche3_state8_seed0.csv")

countMat380 <- as.matrix(omics380)

dim(countMat380)

rownames(countMat380) <- rownames(meta380)


spatial_locs <- data.frame(
  x = meta380$col,
  y = meta380$row,
  row.names = rownames(countMat380)
)

meta380 <- as.data.frame(meta380)

rownames(spatial_locs) <- rownames(countMat380)

countMat380 <- t(countMat380)


giotto_obj <- createGobject(
  countMat = countMat380,
  spatial_locs = spatial_locs,
  cell_metadata = meta380,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "sim/svg_giotto_kmeans380.csv")

fwrite(svg_giotto_rank, file = "sim/svg_giotto_rank380.csv")