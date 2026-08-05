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


meta391 <- read_csv("dataset_0515/sim_meta_niche3_state9_seed1.csv")

omics391 <- read_csv("dataset_0515/sim_omics_niche3_state9_seed1.csv")

countMat391 <- as.matrix(omics391)

dim(countMat391)

rownames(countMat391) <- rownames(meta391)


spatial_locs <- data.frame(
  x = meta391$col,
  y = meta391$row,
  row.names = rownames(countMat391)
)

meta391 <- as.data.frame(meta391)

rownames(spatial_locs) <- rownames(countMat391)

countMat391 <- t(countMat391)


giotto_obj <- createGobject(
  countMat = countMat391,
  spatial_locs = spatial_locs,
  cell_metadata = meta391,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "sim/svg_giotto_kmeans391.csv")

fwrite(svg_giotto_rank, file = "sim/svg_giotto_rank391.csv")