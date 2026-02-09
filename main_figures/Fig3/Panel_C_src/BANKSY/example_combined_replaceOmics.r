library(Banksy)

library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)

library(scater)
library(cowplot)
library(ggplot2)

get_script_dir <- function() {
  # Try to get script directory from commandArgs first (for Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- args[grep(file_arg, args)]
  
  if (length(script_path) > 0) {
    # Running via Rscript
    path <- sub(file_arg, "", script_path)
    return(normalizePath(dirname(path), winslash = "/", mustWork = FALSE))
  } else {
    # Running in RStudio or interactive session
    if (exists("rstudioapi") && rstudioapi::isAvailable()) {
      # Use rstudioapi if available
      return(dirname(rstudioapi::getActiveDocumentContext()$path))
    } else {
      # Fallback: try to load rstudioapi
      if (requireNamespace("rstudioapi", quietly = TRUE)) {
        return(dirname(rstudioapi::getActiveDocumentContext()$path))
      } else {
        # Last resort: use current working directory
        warning("Cannot determine script path. Using current working directory.")
        return(getwd())
      }
    }
  }
}

SCRIPT_DIR <- get_script_dir()

# Load the data
{
  set.seed(0)
  
  real_rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/Xenium_reference_count.csv'),
                            header = T,
                            row.names = 1)
  real_meta <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/Xenium_reference_metadata.csv'),
                        header = T)
  cell_count <- sort(table(real_meta$Cluster), decreasing = T)
  real_meta$state_rank <- 0
  for (i in 1:length(cell_count)) {
    real_meta$state_rank[real_meta$Cluster == names(cell_count)[i]] <- i 
  }
  real_location <- as.data.frame(real_meta[,c("x_centroid", "y_centroid", "Cluster")])
  # normalize the coordinates in real_location so that they sit between 0 and 100 for each column seperately
  real_location$x_centroid <- (real_location[,1] - min(real_location[,1])) / (max(real_location[,1]) - min(real_location[,1])) * 100
  real_location$y_centroid <- (real_location[,2] - min(real_location[,2])) / (max(real_location[,2]) - min(real_location[,2])) * 100
  colnames(real_location) <- c("row", "col", "Cluster")
  real_location$row <- real_location$row + 200
  rownames(real_location) <- paste0("real_", seq(1, nrow(real_location)))
  colnames(real_rawcount) <- rownames(real_location)
  real_rawcount <- as.matrix(real_rawcount)
  
  
  rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/simspace_fitted_count.csv'), 
                       header = T,
                       row.names = 1)
  location <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/simspace_fitted_metadata.csv'),
                       header = T)
  rawcount <- t(rawcount)
  location <- as.matrix(location)
  
  rownames(location) <- colnames(rawcount)
  new_location <- location[,c("row", "col", 'fitted_celltype')]
  new_location <- as.data.frame(new_location)
  ## Match the SimSpace simulated cells with real gene data
  {
    sim_count_new <- matrix(NA, nrow = nrow(real_rawcount), ncol = ncol(real_rawcount))
    rownames(sim_count_new) <- rownames(real_rawcount)
    colnames(sim_count_new) <- colnames(real_rawcount)
    for (cell_type in unique(real_location$Cluster)) {
      # Filter cells with the same state_rank in new_location
      matching_cells <- rownames(real_location)[real_location[, "Cluster"] == cell_type]
      
      # Randomly sample rows from rawcount corresponding to the matching cells
      num_samples <- sum(real_location$Cluster == cell_type)
      sampled_rows <- real_rawcount[, matching_cells, drop = FALSE][, sample(1:length(matching_cells), 
                                                                             size = num_samples, 
                                                                             replace = TRUE), drop = FALSE]
      
      sim_count_new[, rownames(real_location)[real_location$Cluster == cell_type]] <- sampled_rows
    }
    
    sim_count_new <- as.matrix(sim_count_new)
  }
  
  ## scCube data loading
  {
    scCube_rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/scCube_fitted_count.csv'), 
                         header = T,
                         row.names = 1)
    scCube_location <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/scCube_fitted_metadata.csv'),
                         header = T)
    scCube_location$point_x <- (scCube_location$point_x  - min(scCube_location$point_x )) / (max(scCube_location$point_x ) - min(scCube_location$point_x )) * 100
    scCube_location$point_y <- (scCube_location$point_y  - min(scCube_location$point_y )) / (max(scCube_location$point_y ) - min(scCube_location$point_y )) * 100
    
    scCube_location$point_x <- scCube_location$point_x + 100
    scCube_location <- as.matrix(scCube_location)
    
    rownames(scCube_location) <- colnames(scCube_rawcount)
    new_scCube_location <- scCube_location[,c("point_x", "point_y", 'Cell_type')]
    new_scCube_location <- as.data.frame(new_scCube_location)
    colnames(new_scCube_location) <- c("row", "col", 'state_rank')
    
    ## Match the sccube simulated cells with real gene data
    {
      sccube_count_new <- matrix(NA, nrow = nrow(real_rawcount), ncol = ncol(real_rawcount))
      rownames(sccube_count_new) <- rownames(real_rawcount)
      colnames(sccube_count_new) <- colnames(real_rawcount)
      scCube_location <- as.data.frame(scCube_location)
      for (cell_type in unique(scCube_location$Cell_type)) {
        # Filter cells with the same state_rank in new_location
        matching_cells <- rownames(real_location)[real_location[, "Cluster"] == cell_type]
        
        # Randomly sample rows from rawcount corresponding to the matching cells
        num_samples <- sum(real_location$Cluster == cell_type)
        sampled_rows <- real_rawcount[, matching_cells, drop = FALSE][, sample(1:length(matching_cells), 
                                                                               size = num_samples, 
                                                                               replace = TRUE), drop = FALSE]
        
        sccube_count_new[, rownames(real_location)[real_location$Cluster == cell_type]] <- sampled_rows
      }
      
      sccube_count_new <- as.matrix(sccube_count_new)
    }
  }
  
  combined_rawcount <- cbind(rawcount, sccube_count_new, sim_count_new)
  new_location <- as.matrix(new_location)
  real_location <- as.matrix(real_location)
  new_scCube_location <- as.matrix(new_scCube_location)
  combined_location <- rbind(new_location, new_scCube_location, real_location)
  combined_location <- as.data.frame(combined_location)
  combined_location$Dataset <- c(rep("SimSpace", nrow(location)), rep("scCube", nrow(scCube_location)), rep("Xenium Reference", nrow(real_location)))
  
  combined_rawcount <- as.matrix(combined_rawcount)
  combined_location <- as.matrix(combined_location)
}

combined_location <- as.data.frame(combined_location)
combined_location$row <- as.numeric(combined_location$row)
combined_location$col <- as.numeric(combined_location$col)

se <- SpatialExperiment(assay = list(counts = combined_rawcount), 
                        spatialCoords = as.matrix(combined_location[,c(1,2)]))

# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)

lambda <- c(0, 0.2, 0.5, 0.8)
k_geom <- c(15, 30)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)

set.seed(1)
se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 1.2)

se <- Banksy::connectClusters(se)

cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

tmp <- as.data.frame(se@int_colData)

combined_location$PCA1 <- tmp$reducedDims.PCA_M1_lam0.2.PC1
combined_location$PCA2 <- tmp$reducedDims.PCA_M1_lam0.2.PC2
combined_location$UMAP1 <- tmp$reducedDims.UMAP_M1_lam0.2.V1
combined_location$UMAP2 <- tmp$reducedDims.UMAP_M1_lam0.2.V2
combined_location$UMAP1_cluster <- se@colData$clust_M1_lam0.2_k50_res1.2
combined_location$UMAP1_0.5 <- tmp$reducedDims.UMAP_M1_lam0.5.V1
combined_location$UMAP2_0.5 <- tmp$reducedDims.UMAP_M1_lam0.5.V2
combined_location$UMAP1_0.8 <- tmp$reducedDims.UMAP_M1_lam0.8.V1
combined_location$UMAP2_0.8 <- tmp$reducedDims.UMAP_M1_lam0.8.V2


set.seed(0)
combined_location <- combined_location[sample(nrow(combined_location), nrow(combined_location)),]

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = fitted_celltype), size=1.2, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Cell Type") +
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave("Banksy_UMAP_type_replaceOmics.png", 
       width = 5, height = 3, dpi = 300)

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = Dataset), size=1.2, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     name = "Source") +
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave("Banksy_UMAP_dataset_replaceOmics.png", 
       width = 5, height = 3, dpi = 300)

{
  ggplot(data = combined_location) + 
    geom_point(aes(x = UMAP1_0.5, y = UMAP2_0.5, color = fitted_celltype), size=1.2, alpha = 0.6, stroke = 0) + 
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                  "#A65628", "#999999", "#F781BF", "#A65628"),
                       name = "Cell Type") +
    theme_bw() + 
    theme(aspect.ratio = 1)
  ggsave("Banksy_UMAP_type_l0.5_replaceOmics.png", 
         width = 5, height = 3, dpi = 300)
  
  ggplot(data = combined_location) + 
    geom_point(aes(x = UMAP1_0.5, y = UMAP2_0.5, color = Dataset), size=1.2, alpha = 0.6, stroke = 0) + 
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                       name = "Source") +
    theme_bw() + 
    theme(aspect.ratio = 1)
  ggsave("Banksy_UMAP_dataset_l0.5_replaceOmics.png", 
         width = 5, height = 3, dpi = 300)
  
  ggplot(data = combined_location) + 
    geom_point(aes(x = UMAP1_0.8, y = UMAP2_0.8, color = fitted_celltype), size=1.2, alpha = 0.6, stroke = 0) + 
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                  "#A65628", "#999999", "#F781BF", "#A65628"),
                       name = "Cell Type") +
    theme_bw() + 
    theme(aspect.ratio = 1)
  ggsave("Banksy_UMAP_type_l0.8_replaceOmics.png", 
         width = 5, height = 3, dpi = 300)
  
  ggplot(data = combined_location) + 
    geom_point(aes(x = UMAP1_0.8, y = UMAP2_0.8, color = Dataset), size=1.2, alpha = 0.6, stroke = 0) + 
    scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                       name = "Source") +
    theme_bw() + 
    theme(aspect.ratio = 1)
  ggsave("Banksy_UMAP_dataset_l0.8_replaceOmics.png", 
         width = 5, height = 3, dpi = 300)
}








