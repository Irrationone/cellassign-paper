# Merge SCEs

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(pheatmap)
library(Matrix)
library(ggrastr)
library(mclust)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Combine SCEs.")
parser$add_argument('--sces', metavar='FILE', type='character', nargs = '+',
                    help="SCEs to merge")
parser$add_argument('--sample_ids', type='character', nargs = '+',
                    help="Sample IDs to annotate with")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for merged SCE")
args <- parser$parse_args()

sce_files <- unlist(args$sce)
sample_ids <- unlist(args$sample_id)

sces <- lapply(seq_along(sce_files), function(i) {
  sce <- readRDS(sce_files[i]) 
  sce$sample_id <- sample_ids[i]
  return(sce)
})

merge_sces <- function(sces) {
  coldata_union <- Reduce(union, lapply(sces, function(x) {
    colnames(colData(x))
  }))
  
  common_genes <- Reduce(intersect, lapply(sces, function(x) {
    rownames(x)
  }))
  
  sces <- lapply(sces, function(x) {
    missing_cols <- setdiff(coldata_union, colnames(colData(x)))
    rowdat_cols <- setdiff(colnames(rowData(x)), c("mean_counts", "log10_mean_counts",
                                                   "n_cells_by_counts", "pct_dropout_by_counts", "total_counts",
                                                   "log10_total_counts"))
    rowData(x) <- rowData(x)[,rowdat_cols]
    colData(x)[,missing_cols] <- NA
    colData(x) <- colData(x)[,coldata_union]
    x <- x[common_genes,]
    return(x)
  })
  
  sce <- Reduce(cbind, sces)
  
  # Recompute size factors
  qclust <- quickCluster(sce, min.size = 30)
  sce <- computeSumFactors(sce, clusters = qclust)
  
  sce$size_factor <- sizeFactors(sce)
  return(sce)
}

sce_merged <- merge_sces(sces)
sce_merged <- normalize(sce_merged)
sce_merged <- runPCA(sce_merged, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
sce_merged <- runTSNE(sce_merged, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
sce_merged <- runUMAP(sce_merged, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)

message("Saving output ...")

saveRDS(sce_merged, args$outfname)

message("Completed.")

