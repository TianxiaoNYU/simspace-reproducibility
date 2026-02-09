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


meta281 <- read_csv("dataset_0515/sim_meta_niche2_state8_seed1.csv")

omics281 <- read_csv("dataset_0515/sim_omics_niche2_state8_seed1.csv")


countMat281 <- as.matrix(omics281)

dim(countMat281)

rownames(countMat281) <- rownames(meta281)


spatial_locs <- data.frame(
  x = meta281$col,
  y = meta281$row,
  row.names = rownames(countMat281)
)

meta281 <- as.data.frame(meta281)

rownames(spatial_locs) <- rownames(countMat281)

countMat281 <- t(countMat281)


giotto_obj <- createGobject(
  countMat = countMat281,
  spatial_locs = spatial_locs,
  cell_metadata = meta281,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "svg_giotto_kmeans281.csv")

fwrite(svg_giotto_rank, file = "svg_giotto_rank281.csv")