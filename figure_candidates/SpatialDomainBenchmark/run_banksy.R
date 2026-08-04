#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Banksy)
    library(Matrix)
    library(SpatialExperiment)
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(3, 4))) {
    stop("Usage: run_banksy.R DATASET_DIR OUTPUT_DIR SEED [LAMBDA]")
}
dataset_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
seed <- as.integer(args[[3]])
lambda <- if (length(args) == 4) as.numeric(args[[4]]) else 0.5
if (!is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("BANKSY lambda must be between zero and one")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
started <- proc.time()[["elapsed"]]

record <- list(
    method = "banksy",
    dataset = basename(dataset_dir),
    method_seed = seed,
    status = "running"
)

write_record <- function() {
    record$wall_time_seconds <<- proc.time()[["elapsed"]] - started
    record$finished_utc <<- format(Sys.time(), tz = "UTC", usetz = TRUE)
    jsonlite::write_json(
        record,
        path = file.path(output_dir, "run.json"),
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null"
    )
}

tryCatch({
    counts_location_gene <- Matrix::readMM(gzfile(file.path(dataset_dir, "counts.mtx.gz")))
    genes <- read.delim(file.path(dataset_dir, "genes.tsv"), check.names = FALSE)
    coordinates <- read.delim(file.path(dataset_dir, "coordinates.tsv"), check.names = FALSE)
    if (nrow(counts_location_gene) != nrow(coordinates) || ncol(counts_location_gene) != nrow(genes)) {
        stop("BANKSY input dimensions are inconsistent")
    }
    counts <- Matrix::t(counts_location_gene)
    rownames(counts) <- genes$gene_id
    colnames(counts) <- coordinates$location_id
    locs <- as.matrix(coordinates[, c("x", "y")])
    rownames(locs) <- coordinates$location_id
    spe <- SpatialExperiment(
        assays = list(counts = counts),
        spatialCoords = locs
    )
    spe <- scuttle::logNormCounts(spe)
    spe <- Banksy::computeBanksy(
        spe,
        assay_name = "logcounts",
        k_geom = 30,
        compute_agf = FALSE,
        seed = seed,
        verbose = TRUE
    )
    spe <- Banksy::runBanksyPCA(
        spe,
        use_agf = FALSE,
        lambda = lambda,
        npcs = 50,
        seed = seed
    )
    spe <- Banksy::clusterBanksy(
        spe,
        use_agf = FALSE,
        lambda = lambda,
        use_pcs = TRUE,
        npcs = 50,
        algo = "mclust",
        mclust.G = 4,
        seed = seed
    )
    cluster_columns <- Banksy::clusterNames(spe)
    if (length(cluster_columns) < 1) {
        stop("BANKSY did not expose a cluster column")
    }
    cluster_column <- tail(cluster_columns, 1)
    assignments <- data.frame(
        location_id = colnames(spe),
        predicted_domain = as.character(colData(spe)[[cluster_column]])
    )
    write.table(
        assignments,
        file.path(output_dir, "assignments.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    record$status <- "success"
    record$software_version <- paste("BANKSY", as.character(packageVersion("Banksy")))
    record$clustering <- paste0("lambda=", format(lambda, trim = TRUE), ", mclust G=4")
    record$lambda <- lambda
    record$selected_gene_count <- nrow(spe)
    record$k_geom <- 30
    record$morphology_used <- FALSE
    record$cluster_column <- cluster_column
    write_record()
}, error = function(error) {
    record$status <<- "failed"
    record$error_type <<- class(error)[[1]]
    record$error_message <<- conditionMessage(error)
    write_record()
    stop(error)
})
