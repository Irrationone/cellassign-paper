#' Unsupervised clustering of an SCE, subsetting by a particular celltype

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
library(cellassign.utils)
library(Seurat)
library(argparse)

parser <- ArgumentParser(description = "Unsupervised clustering on an SCE")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--celltypes', type = 'character', nargs = '+',
                    help="Name of celltypes to subset on", default = NULL)
parser$add_argument('--clustering_methods', type = 'character', nargs ='+',
                    help="Clustering method(s) to run")
parser$add_argument('--clustering_method_use', type = 'character', nargs ='+',
                    help="Clustering method to use as final annotations")
parser$add_argument('--random_seed', type = 'integer',
                    help="Random seed to use", default = 19279)
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment to use", default = "r-tensorflow")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cluster assignments.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

# Subset on celltype
if (!is.null(args$celltypes)) {
  sce <- sce %>%
    scater::filter(celltype_full %in% args$celltypes)
}

for (method in unlist(args$clustering_methods)) {
  set.seed(args$random_seed)
  if (str_detect(method, "^seurat")) {
    resolution <- str_extract(method, "(?<=_)[0-9\\.]+$")
  } else {
    resolution <- 0.8
  }
  sce <- cluster_wrapper(sce, 
                         gene_subset = NULL, 
                         dimreduce_method = "PCA", 
                         clustering_method = method,
                         seurat_resolution = resolution)
}

cnames <- colnames(colData(sce))
cluster_cols <- cnames[str_detect(cnames, "_cluster$")]
cluster_assignments <- colData(sce)[,c("sample_barcode", cluster_cols)]

# Label final cluster
final_method <- str_replace(args$clustering_method_use, "_[0-9\\.]+$", "")
cluster_assignments[,"cluster"] <- colData(sce)[,paste0(final_method, "_cluster")]

# Write outputs
write.table(cluster_assignments, file = args$outfname, quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE)

cat("Completed.\n")


