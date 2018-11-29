#' Cell type assignments with CellAssign on a SingleCellExperiment

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(microbenchmark)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Run CellAssign on a SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--fc_percentile', type = 'double',
                    help="FC percentile to filter at", default = 0.5)
parser$add_argument('--expr_percentile', type = 'double',
                    help="Expression percentile to filter at", default = 0.5)
parser$add_argument('--max_genes', type = 'integer',
                    help="Maximum number of genes to choose per class", default = 5)
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--seed', type = 'integer',
                    help="Random seed", default = 1483)
parser$add_argument('--max_batch_size', type = 'integer',
                    help="Maximum batch size.", default = 5000)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cellassign.")
parser$add_argument('--outtime', type = 'character', metavar = 'FILE',
                    help="Output path for benchmarking.")
args <- parser$parse_args()

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')
reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

conda_env <- args$conda_env # Tensorflow virtualenv path
fc_percentile <- args$fc_percentile
expr_percentile <- args$expr_percentile
max_genes <- args$max_genes

sce_path <- args$sce
sce <- readRDS(sce_path)

markers_to_use <- select_markers(sce, percentile_fc = fc_percentile, 
                                 percentile_meanexpr = expr_percentile, 
                                 frac_genes = 1, max_genes_per_class = max_genes)

print(markers_to_use)
print(length(markers_to_use))

rho <- create_rho_matrix(sce, markers_to_use, 
                         wrong_proportion = 0)

print(dim(rho))

s <- sizeFactors(sce)
B <- 10
X <- NULL

sce_specific <- sce[rownames(rho),]
counts(sce_specific) <- as.matrix(counts(sce_specific))

set.seed(args$seed)

n_batches <- ceiling(ncol(sce_specific)/args$max_batch_size)

print(n_batches)

timestats <- microbenchmark(res <- cellassign_em(exprs_obj = sce_specific, 
                                                 s = s, 
                                                 rho = rho, 
                                                 X = X, 
                                                 B = B, 
                                                 use_priors = FALSE, 
                                                 prior_type = "shrinkage", 
                                                 delta_variance_prior = FALSE, 
                                                 verbose = FALSE, 
                                                 em_convergence_thres = 1e-3, 
                                                 min_delta = 0, 
                                                 num_runs = 1,
                                                 n_batches = n_batches), times = 1)

# Write outputs
saveRDS(timestats, file = args$outtime)
saveRDS(res, file = args$outfname)

cat("Completed.\n")


