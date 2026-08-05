library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(BiocGenerics)
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
library(scater)

source("call_variable_genes.R")

meta280 <- read_csv("dataset_0515/sim_meta_niche2_state8_seed0.csv")

omics280 <- read_csv("dataset_0515/sim_omics_niche2_state8_seed0.csv")

countMat280 <- as.matrix(omics280)

dim(countMat280)

rownames(countMat280) <- rownames(meta280)

cell_ids <- rownames(countMat280)

#countMat280 <- t(countMat280)

#cell_ids <- colnames(countMat280)

#rownames(countMat280) <- cell_ids
#rownames(meta280) <- cell_ids


spatial_locs <- data.frame(
  x = meta280$col,
  y = meta280$row,
  row.names = rownames(countMat280)
)

meta280 <- as.data.frame(meta280)

rownames(spatial_locs) <- rownames(countMat280)

colnames(countMat280) <- paste0("cell_", 1:ncol(countMat280))
cell_ids <- paste0("", 1:nrow(countMat280))

rownames(spatial_locs) <- cell_ids


countMat280 <- t(countMat280)

# Giotto

giotto_obj <- createGobject(
  countMat = countMat280,
  spatial_locs = spatial_locs,
  cell_metadata = meta280,
  filter = FALSE,  # Whether to perform filtering
  normalize = TRUE  # Whether to perform normalization
)

svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "svg_giotto_kmeans_sim280.csv")

fwrite(svg_giotto_rank, file = "svg_giotto_rank_sim280.csv")