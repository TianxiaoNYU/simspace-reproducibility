library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(viridis)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(path) == 0) stop("Cannot determine script path (not running via Rscript --file=...)")
  normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- get_script_dir()

raw_data <- read.csv(file.path(SCRIPT_DIR, '../Panel_B_C_D_data/fitted_pvals.csv')) 
rownames(raw_data) <- raw_data[,1]
raw_data <- raw_data[, -1]
# Change the values in raw_data into their ranks starting with 1
raw_data <- apply(raw_data, 2, rank, ties.method = "average")
row_means <- rowMeans(raw_data)
# Sort raw_data by the values in the first column
raw_data <- raw_data[order(-row_means), ]
colnames(raw_data) <- c('scBSP', 'Hotspot', 'Giotto(kmeans)', 'Giotto(rank)', 'Moran\'s I', 'SpatialDE', 'SPARK')
pheatmap(raw_data, 
         angle_col=315,
         color = viridis(256),
         # scale = 'row',
         cluster_rows = F,
         cluster_cols = T,
         treeheight_row = 0,
         treeheight_col = 0,
         fontsize = 8,
         show_rownames = F,
         clustering_distance_rows = "euclidean",
         filename = file.path(SCRIPT_DIR, '../Fig5_panel_B.png'),
         width = 2,
         height = 4,
         )

