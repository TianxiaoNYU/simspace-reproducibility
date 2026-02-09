library(tidyverse)
library(SingleCellExperiment)
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


meta291 <- read_csv("dataset_0515/sim_meta_niche2_state9_seed1.csv")

omics291 <- read_csv("dataset_0515/sim_omics_niche2_state9_seed1.csv")

countMat291 <- as.matrix(omics291)

dim(countMat291)

rownames(countMat291) <- rownames(meta291)


spatial_locs <- data.frame(
  x = meta291$col,
  y = meta291$row,
  row.names = rownames(countMat291)
)

meta291 <- as.data.frame(meta291)

rownames(spatial_locs) <- rownames(countMat291)

countMat291 <- t(countMat291)


giotto_obj <- createGobject(
  countMat = countMat291,
  spatial_locs = spatial_locs,
  cell_metadata = meta291,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "sim/svg_giotto_kmeans291.csv")

fwrite(svg_giotto_rank, file = "sim/svg_giotto_rank291.csv")
