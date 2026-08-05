#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
meta_path <- args[1]
omics_path <- args[2]
ref_meta_path <- args[3]
ref_omics_path <- args[4]
save_path <- args[5]
if (length(args) != 5) {
    stop("Usage: Rscript Seurat.R <meta_path> <omics_path> <ref_meta_path> <ref_omics_path> <save_path>")
}

# Load required packages
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
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

reference <- CreateSeuratObject(counts = ref_counts, project = "ref")
reference <- AddMetaData(reference, metadata = cell_types)
reference <- SCTransform(reference, verbose = FALSE) %>%
    RunPCA(verbose = FALSE)

spatial <- CreateSeuratObject(counts = counts, project = "spatial")
spatial <- AddMetaData(spatial, metadata = coords)
spatial <- SCTransform(spatial, verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = reference, query = spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference$Cluster, prediction.assay = TRUE,
    weight.reduction = spatial[["pca"]], dims = 1:30)
res <- predictions.assay@data
res <- as.data.frame(t(res))
res <- res[, -ncol(res)]

res <- res[, sort(colnames(res))]
res$X <- coords$x
res$Y <- coords$y
write.table(res, 
            paste0(save_path, "/Seurat_res.csv"),
            sep = ",",
            quote = F)
