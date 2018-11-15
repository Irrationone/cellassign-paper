# Batch correct data with scanorama/MNN, etc. for visualization

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

parser <- ArgumentParser(description = "Determine malignant status of B cells in FL data")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to follicular SingleCellExperiment RDS")
parser$add_argument('--conda_env', type='character',
                    help="Conda environment with scanorama", default = "r-tensorflow")
parser$add_argument('--method', type='character',
                    help="Batch correction method")
parser$add_argument('--batch_column', type='character',
                    help="Batch column", default = "patient")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for batch corrected SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
batch_col <- args$batch_column

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

# Select features

batches <- unique(colData(sce)[,batch_col])

fit <- trendVar(sce, parametric=TRUE, use.spikes = FALSE)
decomp <- decomposeVar(sce, fit)
decomp$Symbol <- rowData(sce)$Symbol

decomp_chosen <- decomp %>% subset(bio > 0)
chosen <- rownames(decomp_chosen)

chosen_features <- rep(list(chosen), length(batches))

# Apply batch correction

if (args$method == "scanorama") {
  norm_matrices <- lapply(batches, function(bat) {
    sce_sub <- sce[,colData(sce)[,batch_col] == bat]
    norm_data <- t(as.matrix(logcounts(sce_sub)))
    return(norm_data)
  })
  
  indexes <- lapply(batches, function(bat) {
    which(colData(sce)[,batch_col] == bat)
  })
  
  scanorama <- reticulate::import('scanorama')
  
  # Integration and batch correction
  integrated_corrected_data <- scanorama$correct(norm_matrices,
                                                 chosen_features,
                                                 return_dimred = TRUE,
                                                 return_dense = TRUE)
  
  scanorama_mat <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = 100)
  scanorama_mat_expr <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = ncol(integrated_corrected_data[[2]][[1]]))
  for (i in seq_along(indexes)) {
    idx <- indexes[[i]]
    scanorama_mat[idx,] <- integrated_corrected_data[[1]][[i]]
    scanorama_mat_expr[idx,] <- integrated_corrected_data[[2]][[i]]
  }
  reducedDim(sce, "scanorama_int") <- scanorama_mat
  
  sce_sel <- sce
  sce_sel <- sce_sel[integrated_corrected_data[[3]],]
  assay(sce_sel, "scanorama") <- t(scanorama_mat_expr)
  pca_res <- pca2(sce_sel, ntop = 1000, ncomponents = 50, exprs_values = "scanorama")
  reducedDim(sce_sel, "PCA2") <- pca_res$x
  
  sce_sel <- runTSNE(sce_sel, use_dimred = "PCA2", n_dimred = 50)
  sce_sel <- runUMAP(sce_sel, use_dimred = "PCA2", n_dimred = 50)
  
  reducedDim(sce, "scanorama_PCA") <- reducedDim(sce_sel, "PCA2")
  reducedDim(sce, "scanorama_TSNE") <- reducedDim(sce_sel, "TSNE")
  reducedDim(sce, "scanorama_UMAP") <- reducedDim(sce_sel, "UMAP")
  sce@metadata$batchcor_pca_res <- pca_res
  sce@metadata$batchcor_genes <- chosen
} else if (args$method == "MNN") {
  chosen_features_merged <- do.call(union, chosen_features)
  stop("MNN not implemented yet.")
} else {
  stop("Unrecognized method.")
}


# Save malignant status
saveRDS(sce, args$outfname)

cat("Completed.\n")



