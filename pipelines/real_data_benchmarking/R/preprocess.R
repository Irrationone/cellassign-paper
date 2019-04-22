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
parser$add_argument('--umap_neighbors', type='integer',
                    help="Number of nearest neighbours for UMAP", default = 15)
parser$add_argument('--umap_min_dist', type = 'double',
                    help="Minimum distance for UMAP", default = 0.1)
parser$add_argument('--nmads', type='character',
                    help="Number of MADs to filter at.", default = "3")
parser$add_argument('--run_dimred', action='store_true',
                    help="Run dimensionality reduction.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")
args <- parser$parse_args()

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

# Write outputs
saveRDS(sce_normalized, file = args$outfname)

cat("Completed.\n")


