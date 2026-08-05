args <- commandArgs(trailingOnly = TRUE)
meta_path <- args[1]
omics_path <- args[2]
save_path <- args[3]
if (length(args) != 3) {
    stop("Usage: Rscript CARD.R <meta_path> <omics_path> <save_path>")
}

# Load required packages
suppressPackageStartupMessages({
    library(STdeconvolve)
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

corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.01)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = 9)
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

deconProp <- as.data.frame(deconProp)
deconProp <- deconProp[, sort(colnames(deconProp))]
deconProp$X <- coords$x
deconProp$Y <- coords$y
write.table(deconProp, 
            paste0(save_path, "/STdeconvolve_res.csv"),
            sep = ",",
            quote = F)