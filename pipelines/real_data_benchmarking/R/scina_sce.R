#' Cell type assignments with SCINA on a SingleCellExperiment

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
library(SCINA)

parser <- ArgumentParser(description = "Run SCINA on a SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character', nargs = '+',
                    help="Path(s) to SingleCellExperiment RDS")
parser$add_argument('--marker_list', metavar='FILE', type = 'character',
                    help="Path to a marker list yaml", default = NULL)
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--sensitivity_cutoff', type = 'double', 
                    help="Sensitivity cutoff", default = 1)
parser$add_argument('--rm_overlap', type = 'integer', 
                    help="Remove overlapping genes", default = 1, choices = c(0,1))
parser$add_argument('--allow_unknown', type = 'integer', 
                    help="Allow unknown celltypes", default = 1, choices = c(0,1))
parser$add_argument('--batch_correct', type = 'character',
                    help="Batch correct multiple SCEs", default = "none")
parser$add_argument('--celltypes', type = 'character', nargs = '+',
                    help="Celltypes to subset", default = "all")
parser$add_argument('--celltype_col', type = 'character',
                    help="Celltype column", default = "celltype")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cell cycle assignments.")
args <- parser$parse_args()

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

batch_correct2 <- function(sce, batch_col, method = "scanorama") {
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
    
    exprs_bc <- assay(sce_sel, "scanorama")
  } else {
    stop("Other methods not implemented.")
  }
  
  return(list(sce=sce, exprs=exprs_bc))
}

sce_paths <- unlist(args$sce)

if (length(sce_paths) == 1) {
  sce <- readRDS(sce_paths[1])
} else {
  sces <- lapply(sce_paths, function(x) {
    sce <- readRDS(x)
    colData(sce)$path <- x
    return(sce)
  })
  
  coldata_union <- Reduce(union, lapply(sces, function(x) {
    colnames(colData(x))
  }))
  
  common_genes <- Reduce(intersect, lapply(sces, function(x) {
    rownames(x)
  }))
  
  sces <- lapply(seq_along(sces), function(i) {
    x <- sces[[i]]
    missing_cols <- setdiff(coldata_union, colnames(colData(x)))
    rowdat_cols <- setdiff(colnames(rowData(x)), c("mean_counts", "log10_mean_counts",
                                                   "n_cells_by_counts", "pct_dropout_by_counts", "total_counts",
                                                   "log10_total_counts"))
    rowData(x) <- rowData(x)[,rowdat_cols]
    colData(x)[,missing_cols] <- NA
    colData(x) <- colData(x)[,coldata_union]
    x <- x[common_genes,]
    x$batch <- paste0("Batch", i)
    return(x)
  })
  
  sce <- Reduce(cbind, sces)
  
  # Recompute size factors
  qclust <- quickCluster(sce, min.size = 30)
  sce <- computeSumFactors(sce, clusters = qclust)
  
  sce$size_factor <- sizeFactors(sce)
}

if (args$batch_correct != "none") {
  bc_res <- batch_correct2(sce, batch_col = args$batch_correct, method = "scanorama")
}

if (all(unlist(args$celltypes) != "all")) {
  sce$celltype <- as.data.frame(colData(sce))[,args$celltype_col]
  
  sce <- sce %>%
    scater::filter(celltype %in% unlist(args$celltypes))
  
  print(sce)
}

# Process marker gene list
marker_list <- read_yaml(args$marker_list)
marker_list <- lapply(marker_list, function(x) {
  get_ensembl_id(intersect(x, rowData(sce)$Symbol), sce)
})

if (args$batch_correct != "none") {
  message("Using batch corrected values")
  expr_mat <- bc_res$exprs
} else {
  expr_mat <- as.matrix(logcounts(sce))
}

scina_res <- SCINA(expr_mat,
                   marker_list,
                   sensitivity_cutoff = args$sensitivity_cutoff,
                   rm_overlap = args$rm_overlap,
                   allow_unknown = args$allow_unknown)

# Write outputs
saveRDS(scina_res, file = args$outfname)

cat("Completed.\n")


