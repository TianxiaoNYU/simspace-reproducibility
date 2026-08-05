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

combined_location <- read.table(
  file.path(SCRIPT_DIR, 'umap_embedding.csv'),
  header = TRUE, 
  sep = ',',
  row.names = 1)
set.seed(0)
combined_location <- combined_location[sample(nrow(combined_location), nrow(combined_location)),]


ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = Cluster), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", 
                                "#A65628", "#999999", "#F781BF", "#A65628"),
                     name = "Cell Type") +
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave(file.path(SCRIPT_DIR, "SEDR_UMAP_type.png"), 
       width = 5, height = 3, dpi = 300)

ggplot(data = combined_location) + 
  geom_point(aes(x = UMAP1, y = UMAP2, color = Dataset), size=1.8, alpha = 0.6, stroke = 0) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"),
                     name = "Source") +
  theme_bw() + 
  theme(aspect.ratio = 1)
ggsave(file.path(SCRIPT_DIR, "SEDR_UMAP_dataset.png"), 
       width = 5, height = 3, dpi = 300)
