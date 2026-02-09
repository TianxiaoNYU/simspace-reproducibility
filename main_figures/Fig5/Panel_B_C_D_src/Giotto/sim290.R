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

meta290 <- read_csv("dataset_0515/sim_meta_niche2_state9_seed0.csv")

omics290 <- read_csv("dataset_0515/sim_omics_niche2_state9_seed0.csv")

countMat290 <- as.matrix(omics290)

dim(countMat290)

rownames(countMat290) <- rownames(meta290)


spatial_locs <- data.frame(
  x = meta290$col,
  y = meta290$row,
  row.names = rownames(countMat290)
)

meta290 <- as.data.frame(meta290)

rownames(spatial_locs) <- rownames(countMat290)

countMat290 <- t(countMat290)


giotto_obj <- createGobject(
  countMat = countMat290,
  spatial_locs = spatial_locs,
  cell_metadata = meta290,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "svg_giotto_kmeans290.csv")

fwrite(svg_giotto_rank, file = "svg_giotto_rank290.csv")