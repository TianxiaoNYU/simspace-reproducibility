library(ggplot2)
library(scBSP)
library(sparseMatrixStats)
library(Matrix)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(path) == 0) stop("Cannot determine script path (not running via Rscript --file=...)")
  normalizePath(dirname(path), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- get_script_dir()

load_path <- file.path(SCRIPT_DIR, '../benchmark_datasets/')
save_path <- file.path(SCRIPT_DIR)
niche_list <- c(2,3)
state_list <- c(8,9)
seed_list <- c(0,1)

for (i in niche_list) {
  for (j in state_list) {
    for (k in seed_list) {
      # Load the datasets
      meta <- read.table(paste(load_path, 'sim_meta_niche', i, '_state', j, '_seed', k, '.csv', sep = ''),
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
      count <- read.table(paste(load_path, 'sim_omics_niche', i, '_state', j, '_seed', k, '.csv', sep = ''),
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
      count <- t(count)
      count <- as(count, "dgCMatrix")
      location <- meta[, c("row", "col")]
      
      # Run scBSP
      P_values <- scBSP(location, count)
      
      # Save the results
      write.csv(P_values, 
                paste(save_path, 'svg_scBSP_', i, j, k, '.csv', sep = ''), 
                row.names = FALSE,
                quote = FALSE)
    }
  }
}



meta <- read.table(paste0(load_path, '/Xenium_ref/tile_23_meta.csv'),
                   header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- read.table(paste0(load_path, '/Xenium_ref/tile_23_count.csv'),
                    header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- as(as.matrix(count), "dgCMatrix")
location <- meta[, c("x_centroid", "y_centroid")]

# Run scBSP
P_values <- scBSP(location, count)

# Save the results
write.csv(P_values, 
          paste(save_path, 'svg_scBSP_reference.csv', sep = ''), 
          row.names = FALSE,
          quote = FALSE)


## Fitted datasets
meta <- read.table(paste0(load_path, '/Xenium_fitted/tile_23_fitted_meta.csv'),
                   header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- read.table(paste0(load_path, '/Xenium_fitted/tile_23_fitted_count.csv'),
                    header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- as(t(count), "dgCMatrix")
location <- meta[, c("row", "col")]

# Run scBSP
P_values <- scBSP(location, count)

# Save the results
write.csv(P_values, 
          paste(save_path, 'svg_scBSP_fitted.csv', sep = ''), 
          row.names = FALSE,
          quote = FALSE)



