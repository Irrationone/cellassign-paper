#' Preprocess SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Preprocess SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--mito_thres', type='double',
                    help="Mitochondrial threshold")
parser$add_argument('--ribo_thres', type='double',
                    help="Ribosomal threshold")
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--umap_neighbors', type='integer',
                    help="Number of nearest neighbours for UMAP", default = 15)
parser$add_argument('--umap_min_dist', type = 'double',
                    help="Minimum distance for UMAP", default = 0.1)
parser$add_argument('--nmads', type='character',
                    help="Number of MADs to filter at.", default = "3")
parser$add_argument('--run_dimred', action='store_true',
                    help="Run dimensionality reduction.")
parser$add_argument('--batch_correct', type = 'character',
                    help="Batch correction variable", default = 'none')
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
args <- parser$parse_args()

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

batch_correct <- function(sce, batch_col, method = "scanorama") {
  ## Feature selection
  batches <- unique(colData(sce)[,batch_col])
  
  fit <- trendVar(sce, parametric=TRUE, use.spikes = FALSE)
  decomp <- decomposeVar(sce, fit)
  decomp$Symbol <- rowData(sce)$Symbol
  
  decomp_chosen <- decomp %>% subset(bio > 0)
  chosen <- rownames(decomp_chosen)
  
  chosen_features <- rep(list(chosen), length(batches))
  
  ## Batch correction
  if (method == "scanorama") {
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
    
    sce <- sce[integrated_corrected_data[[3]],]
    
    print(dim(sce))
    print(dim(assay(sce_sel, "scanorama")))
    
    logcounts(sce) <- assay(sce_sel, "scanorama")
  } else {
    stop("Other methods not implemented.")
  }
  
  return(sce)
}


sce_path <- args$sce
sce <- readRDS(sce_path)

nmads <- as.numeric(args$nmads)

if (args$mito_thres >= 100 & args$ribo_thres >= 100 & is.infinite(nmads)) {
  ## Don't preprocess for dummy thresholds
  sce_filtered <- sce
} else {
  # Get ensembl gene IDs of mitochondrial and ribosomal genes
  mito_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^MT\\-") & !is.na(rowData(sce)$Symbol)]) %>% 
    get_ensembl_id(sce)
  
  ribo_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)") & !is.na(rowData(sce)$Symbol)]) %>%
    get_ensembl_id(sce)
  
  sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
                              list(mitochondrial=mito_genes, ribosomal=ribo_genes))
  
  sce_filtered <- filter_cells(sce, nmads = nmads, type = "lower", 
                               log = TRUE, max_mito = args$mito_thres, max_ribo = args$ribo_thres)
}

if (is.null(sizeFactors(sce_filtered))) {
  qclust <- quickCluster(sce_filtered, min.size = 30)
  sce_filtered <- computeSumFactors(sce_filtered, clusters = qclust)
}

sce_filtered$size_factor <- sizeFactors(sce_filtered)

sce_normalized <- normalize(sce_filtered)

if (args$run_dimred) {
  sce_normalized <- runPCA(sce_normalized, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
  sce_normalized <- runTSNE(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
  
  umap_config <- umap:::umap.defaults
  umap_config$n_neighbors <- args$umap_neighbors
  umap_config$min_dist <- args$umap_min_dist
  
  sce_normalized <- runUMAP(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2, config = umap_config)
}

if (args$batch_correct != "none") {
  sce <- batch_correct(sce, batch_col = args$batch_correct, method = "scanorama")
}

# Write outputs
saveRDS(sce_normalized, file = args$outfname)

cat("Completed.\n")


