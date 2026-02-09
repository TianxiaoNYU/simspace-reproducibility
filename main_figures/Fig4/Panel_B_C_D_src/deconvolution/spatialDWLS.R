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
    library(Giotto)
    library(TOAST)
})


counts <- read.table(omics_path,
                     sep = ",",
                     header = T,
                     row.names = NULL)
counts <- t(counts)
colnames(counts) <- paste0("spot_", 1:ncol(counts))
coords_raw <- read.table(
    meta_path,
    sep = ",",
    header = T,
    row.names = NULL)
coords <- coords_raw[,1:2]
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

# Compute mean expression per cell type
cell_type_labels <- cell_types$Cluster
unique_types <- unique(cell_type_labels)
mean_expr_mat <- sapply(unique_types, function(ct) {
    cells <- which(cell_type_labels == ct)
    if (length(cells) == 1) {
        ref_counts[, cells]
    } else {
        rowMeans(ref_counts[, cells, drop = FALSE])
    }
})
colnames(mean_expr_mat) <- unique_types



instrs <- createGiottoInstructions(python_path = '/Users/zhaotianxiao/mambaforge/envs/simulation/bin/python',
                                   show_plot = F, return_plot = T, save_plot = T,
                                   dpi = 300, height = 9, width = 9)
visium_brain <- createGiottoObject(raw_exprs = counts,
                                   spatial_locs = coords,
                                   instructions = instrs,
                                   cell_metadata = coords_raw)
visium_brain <- filterGiotto(gobject = visium_brain,
                             expression_threshold = 1,
                             gene_det_in_min_cells = 1,
                             min_det_genes_per_cell = 1,
                             expression_values = c('raw'),
                             verbose = T)
visium_brain <- normalizeGiotto(gobject = visium_brain)
gene_metadata <- fDataDT(visium_brain)
visium_brain <- runPCA(gobject = visium_brain, scale_unit = F)
visium_brain <- runUMAP(visium_brain, dimensions_to_use = 1:10)
visium_brain <- createNearestNetwork(gobject = visium_brain, dimensions_to_use = 1:10, k = 15)
visium_brain <- doLeidenCluster(gobject = visium_brain, resolution = 0.4, n_iterations = 1000)


visium_brain <- runDWLSDeconv(visium_brain, sign_matrix = mean_expr_mat, n_cell = 20)
deconv_res <- as.data.frame(visium_brain@spatial_enrichment$DWLS)
deconv_res <- deconv_res[,-1]

deconv_res <- deconv_res[, sort(colnames(deconv_res))]
deconv_res$X <- coords$x
deconv_res$Y <- coords$y
write.table(deconv_res, 
            paste0(save_path, "/spatialDWLS_res.csv"),
            sep = ",",
            quote = F)