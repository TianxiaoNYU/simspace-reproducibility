#!/usr/bin/env Rscript

# Generate gene-level summaries for ten independently seeded molecular
# realizations using one scDesign3 fit to the Figure 4 Xenium tile.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scDesign3)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript generate_molecular_replicates.R",
      "<reference_metadata.csv> <reference_counts.csv>",
      "<molecular_simulation_design.tsv> <output.tsv>"
    )
  )
}

reference_metadata_path <- args[1]
reference_counts_path <- args[2]
simulation_design_path <- args[3]
output_path <- args[4]

reference_metadata <- read.csv(
  reference_metadata_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
reference_counts <- read.csv(
  reference_counts_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
simulation_design <- read.delim(
  simulation_design_path,
  header = TRUE,
  check.names = FALSE
)

if (nrow(reference_metadata) != ncol(reference_counts)) {
  stop("Reference metadata rows and count-matrix columns do not align.")
}
required_design_columns <- c(
  "molecular_seed", "cell_index", "row", "col", "fitted_celltype"
)
if (!all(required_design_columns %in% colnames(simulation_design))) {
  stop("Molecular simulation design is missing required columns.")
}

reference_metadata$celltype <- reference_metadata$Cluster
reference_metadata$row <- reference_metadata$x_centroid
reference_metadata$col <- reference_metadata$y_centroid

sce <- SingleCellExperiment(
  list(counts = as.matrix(reference_counts)),
  colData = reference_metadata
)
logcounts(sce) <- log1p(counts(sce))

set.seed(0)
message("Fitting one scDesign3 molecular model to the Xenium reference...")
model_data <- construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "celltype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("row", "col"),
  corr_by = "celltype"
)
marginal_fit <- fit_marginal(
  data = model_data,
  predictor = "gene",
  mu_formula = "celltype",
  sigma_formula = "celltype",
  family_use = "nb",
  n_cores = 6,
  usebam = FALSE,
  parallelization = "pbmcmapply"
)
copula_fit <- fit_copula(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal_fit,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 6,
  input_data = model_data$dat
)

summaries <- vector("list", 10)
for (seed in 0:9) {
  message(sprintf("Generating molecular realization for seed %d...", seed))
  new_metadata <- simulation_design[
    simulation_design$molecular_seed == seed,
    c("row", "col", "fitted_celltype"),
    drop = FALSE
  ]
  colnames(new_metadata)[3] <- "celltype"
  new_metadata$corr_group <- new_metadata$celltype

  parameters <- extract_para(
    sce = sce,
    marginal_list = marginal_fit,
    n_cores = 6,
    family_use = "nb",
    new_covariate = new_metadata,
    data = model_data$dat
  )
  set.seed(seed)
  simulated_counts <- simu_new(
    sce = sce,
    mean_mat = parameters$mean_mat,
    sigma_mat = parameters$sigma_mat,
    zero_mat = parameters$zero_mat,
    quantile_mat = NULL,
    copula_list = copula_fit$copula_list,
    n_cores = 6,
    family_use = "nb",
    input_data = model_data$dat,
    new_covariate = new_metadata,
    important_feature = copula_fit$important_feature,
    filtered_gene = model_data$filtered_gene
  )
  simulated_counts <- as.matrix(simulated_counts)
  if (nrow(simulated_counts) != nrow(reference_counts)) {
    stop(sprintf("Unexpected gene count for molecular seed %d.", seed))
  }
  gene_names <- rownames(simulated_counts)
  if (is.null(gene_names)) {
    gene_names <- rownames(reference_counts)
  }
  summaries[[seed + 1]] <- data.frame(
    molecular_seed = seed,
    gene = gene_names,
    mean_count = rowMeans(simulated_counts),
    variance = apply(simulated_counts, 1, var),
    zero_fraction = rowMeans(simulated_counts == 0),
    n_cells = ncol(simulated_counts),
    check.names = FALSE
  )
}

summary_table <- do.call(rbind, summaries)
write.table(
  summary_table,
  output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
software_versions <- data.frame(
  software = c("R", "scDesign3", "SingleCellExperiment"),
  version = c(
    as.character(getRversion()),
    as.character(packageVersion("scDesign3")),
    as.character(packageVersion("SingleCellExperiment"))
  )
)
write.table(
  software_versions,
  file.path(dirname(output_path), "molecular_software_versions.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message(sprintf("Wrote %s", output_path))
