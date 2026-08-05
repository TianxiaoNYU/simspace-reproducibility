library(ggplot2)
library(SPARK)

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

rep_list <- c(280, 281, 290, 291, 380, 381, 390, 391)

for (i in niche_list) {
  for (j in state_list) {
    for (k in seed_list) {
      if (i == 2) {
        if (j == 9 && k ==1) {
        }else{
          next
        }
      }
      # Load the datasets
      print(paste("Processing niche:", i, "state:", j, "seed:", k))
      meta <- read.table(paste(load_path, 'sim_meta_niche', i, '_state', j, '_seed', k, '.csv', sep = ''),
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
      count <- read.table(paste(load_path, 'sim_omics_niche', i, '_state', j, '_seed', k, '.csv', sep = ''),
                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
      count <- as.data.frame(t(count))
      location <- meta[, c("row", "col")]
      
      rownames(location) <- colnames(count)
      spark <- CreateSPARKObject(counts=count, 
                                 location=location,
                                 percentage = 0, 
                                 min_total_counts = 10)
      spark@lib_size <- apply(spark@counts, 2, sum)
      spark <- spark.vc(spark, 
                        covariates = NULL, 
                        lib_size = spark@lib_size, 
                        num_core = 10,
                        verbose = F)
      spark <- spark.test(spark, 
                          check_positive = T, 
                          verbose = F)
      res <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
      # Save the results
      write.csv(res, 
                paste(save_path, 'svg_spark_', i, j, k, '.csv', sep = ''), 
                row.names = TRUE,
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
rownames(location) <- colnames(count)
count <- as.matrix(count)

spark <- CreateSPARKObject(counts=count, 
                           location=location,
                           percentage = 0, 
                           min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 8,
                  verbose = F)
spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = F)
res <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]

# Save the results
write.csv(res, 
          paste(save_path, 'svg_spark_reference.csv', sep = ''), 
          row.names = TRUE,
          quote = FALSE)


## Fitted datasets
meta <- read.table(paste0(load_path, '/Xenium_fitted/tile_23_fitted_meta.csv'),
                   header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- read.table(paste0(load_path, '/Xenium_fitted/tile_23_fitted_count.csv'),
                    header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
count <- as(t(count), "dgCMatrix")
location <- meta[, c("row", "col")]

rownames(location) <- colnames(count)
spark <- CreateSPARKObject(counts=count, 
                           location=location,
                           percentage = 0, 
                           min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 8,
                  verbose = F)
spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = F)
res <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]

# Save the results
write.csv(res, 
          paste(save_path, 'svg_spark_fitted.csv', sep = ''), 
          row.names = TRUE,
          quote = FALSE)