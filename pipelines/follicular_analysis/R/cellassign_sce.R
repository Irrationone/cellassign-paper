#' Cell type assignments with CellAssign on a SingleCellExperiment

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Run CellAssign on a SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--marker_gene_matrix', metavar='FILE', type = 'character',
                    help="Path to a rho matrix", default = NULL)
parser$add_argument('--marker_list', metavar='FILE', type = 'character',
                    help="Path to a marker list yaml", default = NULL)
parser$add_argument('--include_other', type = 'character',
                    help="Include an other column in the marker gene matrix if using a marker_list")
parser$add_argument('--num_runs', type = 'double',
                    help="Number of runs", default = 1)
parser$add_argument('--rbf_pieces', type = 'integer',
                    help="Number of pieces to fit (RBF)", default = 20)
parser$add_argument('--min_delta', type = 'double',
                    help="Minimum delta value", default = log(2))
parser$add_argument('--delta_prior', action = 'store_true',
                    help="Use delta shrinkage prior")
parser$add_argument('--variance_prior', action = 'store_true',
                    help="Use delta variance prior")
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cell cycle assignments.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

# Process marker gene matrix
if (is.null(args$marker_gene_matrix)) {
  if (is.null(args$marker_list)) {
    stop("Must specify either a marker gene matrix or a marker list file.")
  }
  
  marker_list <- read_yaml(args$marker_list)
  rho <- marker_list_to_mat(marker_list, include_other = as.logical(args$include_other))
} else {
  rho <- read.table(args$marker_gene_matrix, 
                    row.names = 1, 
                    header = TRUE, 
                    sep = "\t",
                    stringsAsFactors = FALSE)
}

rho <- as.matrix(rho)
s <- sizeFactors(sce)

sce_markers <- sce[get_ensembl_id(rownames(rho), sce),]
rownames(sce_markers) <- rownames(rho)
counts(sce_markers) <- as.matrix(counts(sce_markers))

cellassign_res <- cellassign_em(exprs_obj = sce_markers, 
                                s = s, 
                                rho = rho, 
                                X = NULL, # no batch effect design matrix
                                B = args$rbf_pieces, 
                                use_priors = args$delta_prior, 
                                prior_type = "shrinkage", 
                                delta_variance_prior = args$variance_prior, 
                                verbose = FALSE, 
                                em_convergence_thres = 1e-5, 
                                num_runs = args$num_runs, 
                                min_delta = args$min_delta)

# Write outputs
saveRDS(cellassign_res, file = args$outfname)

cat("Completed.\n")


