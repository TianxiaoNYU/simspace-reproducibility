args <- commandArgs(trailingOnly = TRUE)
nrow <- as.numeric(args[1])
ncol <- as.numeric(args[2])
seed <- as.numeric(args[3])
num_genes <- if (length(args) >= 4) as.integer(args[4]) else 220L

if (is.na(num_genes) || num_genes <= 0) {
  stop("num_genes must be a positive integer")
}

library(scMultiSim)
save_path <- "./"

lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(5.2, 5.9)
)

data(GRN_params_100)
set.seed(seed)

options_ <- list(
  GRN = GRN_params_100,
  speed.up = TRUE,
  num.genes = num_genes,
  num.cells = nrow * ncol,
  num.cifs = 20,
  cif.sigma = 0.2,
  tree = Phyla3(),
  intrinsic.noise = 0.5,
  cci = list(
    params = lig_params,
    max.neighbors = 4,
    grid.size = 30,
    cell.type.interaction = "random",
    step.size = 0.5
  )
)
results <- sim_true_counts(options_)












