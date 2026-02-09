#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
meta_path <- args[1]
omics_path <- args[2]
ref_meta_path <- args[3]
ref_omics_path <- args[4]
save_path <- args[5]
if (length(args) != 5) {
    stop("Usage: Rscript CARD.R <meta_path> <omics_path> <ref_meta_path> <ref_omics_path> <save_path>")
}

# Load required packages
suppressPackageStartupMessages({
    library(CARD)
    library(TOAST)
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
colnames(counts) <- rownames(coords)

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
rownames(cell_types) <- colnames(ref_counts)
cell_types$Sample <- 1

CARD_obj = createCARDObject(
  sc_count = ref_counts,
  sc_meta = cell_types,
  spatial_count = counts,
  spatial_location = coords,
  ct.varname = "Cluster",
  ct.select = sort(unique(cell_types$Cluster)),
  sample.varname = "Sample",
  minCountGene = 1,
  minCountSpot = 1) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

CARD_res <- CARD_obj@Proportion_CARD
CARD_res <- as.data.frame(CARD_res)
CARD_res <- CARD_res[, sort(colnames(CARD_res))]
CARD_res$X <- coords$x
CARD_res$Y <- coords$y
write.table(CARD_res, 
            paste0(save_path, "/CARD_res.csv"),
            sep = ",",
            quote = F)