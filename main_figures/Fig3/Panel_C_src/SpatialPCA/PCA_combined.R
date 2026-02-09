library(SpatialPCA)
library(ggplot2)
library(tidyr)
library(dplyr)

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

{
  rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/simspace_fitted_count.csv'), 
                       header = T,
                       row.names = 1)
  location <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/simspace_fitted_metadata.csv'),
                       header = T,
                       row.names = 1)
  rawcount <- t(rawcount)
  location <- as.matrix(location)
  
  rownames(location) <- colnames(rawcount)
  new_location <- location[,c("row", "col")]
  
  real_rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/Xenium_reference_count.csv'),
                            header = T,
                            row.names = 1)
  real_meta <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/Xenium_reference_metadata.csv'),
                        header = T)
  real_location <- (real_meta[,c("x_centroid", "y_centroid")])
  # normalize the coordinates in real_location so that they sit between 0 and 100 for each column seperately
  real_location$x_centroid <- (real_location[,1] - min(real_location[,1])) / (max(real_location[,1]) - min(real_location[,1])) * 100
  real_location$y_centroid <- (real_location[,2] - min(real_location[,2])) / (max(real_location[,2]) - min(real_location[,2])) * 100
  colnames(real_location) <- c("row", "col")
  real_location$row <- real_location$row + 200
  rownames(real_location) <- paste0("real_", seq(1, nrow(real_location)))
  colnames(real_rawcount) <- rownames(real_location)
  real_rawcount <- as.matrix(real_rawcount)
  
  ## Load scCube data
  {
    scCube_rawcount <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/scCube_fitted_count.csv'), 
                         header = T,
                         row.names = 1)
    scCube_location <- read.csv(file.path(SCRIPT_DIR, '../../Panel_B_C_D_data/scCube_fitted_metadata.csv'),
                         header = T,
                         row.names = 1)
    scCube_location$point_x <- scCube_location$point_x + 100
    scCube_location <- as.matrix(scCube_location)
    
    rownames(scCube_location) <- colnames(scCube_rawcount)
    
    new_scCube_location <- scCube_location[,c("point_x", "point_y")]
    colnames(new_scCube_location) <- c("row", "col")
  }
  
  combined_rawcount <- cbind(rawcount, real_rawcount, scCube_rawcount)
  combined_location <- rbind(new_location, real_location, new_scCube_location)
  
  combined_rawcount <- as.matrix(combined_rawcount)
  combined_location <- as.matrix(combined_location)
  class(combined_location) <- "numeric"
}

# ST = CreateSpatialPCAObject(
#   counts=combined_rawcount, 
#   location=combined_location, 
#   project = "SpatialPCA",
#   gene.type="spatial",
#   sparkversion="spark", 
#   gene.number=100,
#   numCores_spark = 8,
#   customGenelist=NULL,
#   min.loctions=1, 
#   min.features=20)
# ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
# ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
# ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
# saveRDS(ST, file.path(SCRIPT_DIR, 'SpatialPCAres.rds'))

ST <- readRDS(file.path(SCRIPT_DIR, 'SpatialPCAres.rds'))
combined_location <- as.data.frame(combined_location)
combined_location$Cluster <- c(location[,5], real_meta$Cluster, scCube_location[,2])
combined_location$Dataset <- c(rep("SimSpace", nrow(location)), rep("Xenium Reference", nrow(real_location)), rep("scCube", nrow(scCube_location)))
ST_PC <- as.data.frame(t(ST@SpatialPCs))
ST_PC$row <- ST@location[,1]
ST_PC$col <- ST@location[,2]
ST_PC <- left_join(ST_PC, combined_location, by = c("row", "col"))

set.seed(0)
ST_PC <- ST_PC[sample(nrow(ST_PC), nrow(ST_PC)),]

ggplot(data = ST_PC) + 
  geom_point(aes(x = V1, y = V2, color = Cluster), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Cell Type") +
  xlab("PC1") +
  ylab("PC2") + 
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave(file.path(SCRIPT_DIR, "PCA_Celltype.png"), 
       width = 5, height = 3, dpi = 300)


ggplot(data = ST_PC) + 
  geom_point(aes(x = V1, y = V2, color = Dataset), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                  "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Source") +
  xlab("PC1") +
  ylab("PC2") + 
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave(file.path(SCRIPT_DIR, "PCA_Dataset.png"), 
       width = 5, height = 3, dpi = 300)



