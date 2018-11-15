# Train scvis

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(scvis)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Train scvis")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to follicular SingleCellExperiment RDS")
parser$add_argument('--conda_env', type='character',
                    help="Conda environment with scanorama", default = "r-tensorflow")
parser$add_argument('--scvis_config_path', metavar='FILE', type='character',
                    help="Path to scvis config")
parser$add_argument('--scvis_dimname', type='character',
                    help="Name of reduced dim to use for scvis")
parser$add_argument('--scvis_output_dir', metavar='DIR', type='character',
                    help="Path to scvis temporary output directory")
parser$add_argument('--outfname', metavar='FILE', type='character',
                    help="Path to output SCE annotated with scvis results")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
scvis_config_path <- args$scvis_config_path
scvis_dimname <- args$scvis_dimname

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

# Rerun PCA
if (scvis_dimname == "PCA") {
  pca_res <- pca2(sce, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
  reducedDim(sce, "PCA2") <- pca_res$x
  
  sce@metadata$pca_res <- pca_res
  scvis_dimname <- "PCA2"
} else {
  sce@metadata$pca_res <- sce@metadata$batchcor_pca_res
}

scvis_trained <- scvis_train(sce, 
                             output_dir = args$scvis_output_dir,
                             config_file = scvis_config_path, 
                             use_reducedDim = TRUE, 
                             reducedDim_name = scvis_dimname)

scvis_trained$scvis_ll <- reducedDim(scvis_trained, "scvis")[,3]
reducedDim(scvis_trained, "scvis") <- reducedDim(scvis_trained, "scvis")[,1:2]

# Save output
saveRDS(scvis_trained, args$outfname)

cat("Completed.\n")



