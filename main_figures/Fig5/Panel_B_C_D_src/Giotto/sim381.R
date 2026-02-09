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


meta381 <- read_csv("dataset_0515/sim_meta_niche3_state8_seed1.csv")

omics381 <- read_csv("dataset_0515/sim_omics_niche3_state8_seed1.csv")

countMat381 <- as.matrix(omics381)

dim(countMat381)

rownames(countMat381) <- rownames(meta381)


spatial_locs <- data.frame(
  x = meta381$col,
  y = meta381$row,
  row.names = rownames(countMat381)
)

meta381 <- as.data.frame(meta381)

rownames(spatial_locs) <- rownames(countMat381)

countMat381 <- t(countMat381)


giotto_obj <- createGobject(
  countMat = countMat381,
  spatial_locs = spatial_locs,
  cell_metadata = meta381,
  filter = FALSE,
  normalize = TRUE
)


svg_giotto_kmeans = callSVG.Giotto(giotto_obj, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(giotto_obj, svg_method = "rank", runHVG = TRUE, n_cores=32)

fwrite(svg_giotto_kmeans, file = "sim/svg_giotto_kmeans381.csv")

fwrite(svg_giotto_rank, file = "sim/svg_giotto_rank381.csv")