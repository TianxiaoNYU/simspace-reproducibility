#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
meta_path <- args[1]
omics_path <- args[2]
ref_meta_path <- args[3]
ref_omics_path <- args[4]
thread_num <- as.numeric(args[5])
save_path <- args[6]

if (length(args) != 6) {
    stop("Usage: Rscript CARD.R <meta_path> <omics_path> <ref_meta_path> <ref_omics_path> <thread_num> <save_path>")
}

# Load required packages
suppressPackageStartupMessages({
    library(spacexr)
    library(Matrix)
    library(doParallel)
})


counts <- read.table(omics_path,
                     sep = ",",
                     header = T,
                     row.names = NULL)
counts <- t(counts)
colnames(counts) <- paste0("spot_", 1:ncol(counts))
coords <- read.table(meta_path,
                     sep = ",",
                     header = T,
                     row.names = NULL)
coords <- coords[,1:2]
colnames(coords) <- c('x', 'y')
rownames(coords) <- colnames(counts)

ref_counts <- read.table(ref_omics_path,
                         sep = ",",
                         header = T,
                         row.names = 1)
# ref_counts <- t(ref_counts)
ref_counts <- as.matrix(ref_counts) 
colnames(ref_counts) <- paste0("cell_", 1:ncol(ref_counts))
cell_types <- read.table(ref_meta_path,
                         sep = ",",
                         header = T,
                         row.names = 1)
cell_types <- cell_types$Cluster
cell_types <- as.character(cell_types) # convert to factor data type
cell_types <- as.factor(cell_types) 
names(cell_types) <- colnames(ref_counts)

reference <- Reference(ref_counts, cell_types, , min_UMI = 5)
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
puck <- SpatialRNA(coords, counts, nUMI)

## Run RCTD in multi mode
myRCTD_multi <- create.RCTD(puck, reference, max_cores = thread_num, CELL_MIN_INSTANCE = 1, UMI_min = min(nUMI), UMI_min_sigma = min(nUMI))
myRCTD_multi <- run.RCTD(myRCTD_multi, doublet_mode = 'multi')
results <- myRCTD_multi@results
weights <- data.frame(matrix(0, length(results), length(levels(cell_types))))
rownames(weights) <- colnames(counts)
colnames(weights) <- levels(cell_types)
for(i in 1:length(results)){
  weight <- results[[i]]$sub_weights
  weights[i, names(weight)] <- weight
}
weights <- as.data.frame(weights)
weights <- weights[, sort(colnames(weights))]
weights$X <- coords$x
weights$Y <- coords$y
write.table(weights, 
            paste0(save_path, "/RCTD_multi_res.csv"),
            sep = ",",
            quote = F)

## Run RCTD in full mode
myRCTD_full <- create.RCTD(puck, reference, max_cores = thread_num, CELL_MIN_INSTANCE = 1, UMI_min = min(nUMI), UMI_min_sigma = min(nUMI))
myRCTD_full <- run.RCTD(myRCTD_full, doublet_mode = 'full')
results <- myRCTD_full@results$weights
results <- as.matrix(results)
rownames(results) <- colnames(counts)

results <- as.data.frame(results)
results <- results[, sort(colnames(results))]
results$X <- coords$x
results$Y <- coords$y
write.table(results, 
            paste0(save_path, "/RCTD_full_res.csv"),
            sep = ",",
            quote = F)