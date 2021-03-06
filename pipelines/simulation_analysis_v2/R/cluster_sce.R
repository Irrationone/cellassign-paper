
library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(SC3)
library(Matrix)
library(mclust)
library(multiClust)
library(limma)
library(edgeR)
library(infotheo)
library(future)
library(tensorflow)
library(crch)
library(splitstackshape)
library(SCINA)

library(scrna.utils)
library(cellassign)
library(cellassign.utils)
library(scrna.sceutils)
library(Seurat)

library(methods)
library(argparse)

parser <- ArgumentParser(description = "Cluster SCE.")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Input SCE object")
parser$add_argument('--outfname', metavar='FILE', type = 'character',
                    help="Output directory for clustering results.")
parser$add_argument('--dimreduce_method', type = 'character',
                    help="Dimensionality reduction method.")
parser$add_argument('--conda_env', type = 'character',
                    help="Conda environment", default = 'r-tensorflow')
parser$add_argument('--fc_percentile', type = 'double',
                    help="Fold change percentile cutoff", default = 0.95)
parser$add_argument('--expr_percentile', type = 'double',
                    help="Mean expression percentile cutoff", default = 0.9)
parser$add_argument('--frac_genes', type = 'double',
                    help="Fraction of genes to use", default = 1)
parser$add_argument('--test_proportion', type = 'double',
                    help="Test proportion (for supervised methods)", default = 0.5)
parser$add_argument('--clust_method', type = 'character',
                    help="Clustering method")
parser$add_argument('--max_genes', type = 'integer',
                    help="Maximum number of marker genes per cell type", default = 15)
parser$add_argument('--marker_setting', type = 'character',
                    help="Full or marker genes?", default = "full")
args <- parser$parse_args()

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

## Helper functions
cluster_wrapper <- function(sce, gene_subset = NULL, dimreduce_method, clustering_method, conda_env = "r-tensorflow", conda_binary = "/home/rstudio/miniconda/bin/conda", seurat_resolution = 0.6,
                            object2 = NULL, object2_cluster_label = "Group") {
  sce_original <- sce
  if (!is.null(gene_subset)) {
    sce <- sce[gene_subset,]
    sce <- runPCA(sce, ntop = 500, ncomponents = 50, exprs_values = "logcounts")
    
    if (!is.null(object2)) {
      object2 <- object2[gene_subset,]
    }
  }
  
  if (clustering_method == "phenograph") {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, conda_env = conda_env, 
                         conda_binary = conda_binary,
                         phenograph_module = "scrnatools.methods.clustering.phenograph_analysis")
  } else if (clustering_method == "dbscan") {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, dbscan_epsilon = 0.5, dbscan_minPoints = 5)
  } else if (str_detect(clustering_method, "^seurat")) {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method,
                         seurat_resolution = seurat_resolution)
  } else if (clustering_method %in% c("Zheng_cor", "scmap")) {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method,
                         object2 = object2, object2_cluster_label = object2_cluster_label)
  } else {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, ncores = 1)
  }
  
  sce_original$cluster <- sce$cluster
  return(sce_original)
}

dimreduce_method <- args$dimreduce_method
fc_percentile <- args$fc_percentile
expr_percentile <- args$expr_percentile
frac_genes <- args$frac_genes
test_proportion <- args$test_proportion
clustering_method <- args$clust_method
max_genes <- args$max_genes
marker_setting <- args$marker_setting

sce_total <- readRDS(args$sce)

## Split SCE into training and test sets
coldat <- as.data.frame(colData(sce_total))
sets <- splitstackshape::stratified(coldat, group = "Group", size = (1-test_proportion), bothSets = TRUE)

rowData(sce_total)$symbol <- rowData(sce_total)$Gene
rowData(sce_total)$feature_symbol <- rowData(sce_total)$Gene
rowData(sce_total)$Symbol <- rowData(sce_total)$Gene

counts(sce_total) <- as.matrix(counts(sce_total))

sce_train <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP1$Cell))]
sce <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP2$Cell))]


## Clustering code
if (marker_setting == "markers") {
  markers_to_use <- select_markers(sce_total, 
                                   percentile_fc = fc_percentile, 
                                   percentile_meanexpr = expr_percentile, 
                                   frac_genes = frac_genes, 
                                   max_genes_per_class = max_genes)
  gene_subset <- markers_to_use
  n_markers <- length(markers_to_use)
} else if (marker_setting == "full") {
  gene_subset <- NULL
  n_markers <- NA
} else {
  stop("Unrecognized marker setting.")
}

if (str_detect(clustering_method, "^seurat_")) {
  seurat_resolution <- as.numeric(str_replace(clustering_method, "^seurat_", ""))

  sce_sim <- cluster_wrapper(sce, gene_subset = gene_subset, dimreduce_method = dimreduce_method, clustering_method = "seurat", seurat_resolution = seurat_resolution)
  
  cluster_df <- data.frame(cell=as.character(colData(sce_sim)$Cell), cluster=sce_sim$cluster)
} else if (clustering_method %in% c("Zheng_cor", "scmap")) {
  sce_sim <- cluster_wrapper(sce, gene_subset = gene_subset, dimreduce_method = dimreduce_method, clustering_method = clustering_method, object2 = sce_train, object2_cluster_label = "Group")
  
  cluster_df <- data.frame(cell=as.character(colData(sce_sim)$Cell), cluster=sce_sim$cluster)
} else {
  sce_sim <- cluster_wrapper(sce, gene_subset = gene_subset, dimreduce_method = dimreduce_method, clustering_method = clustering_method)

  cluster_df <- data.frame(cell=as.character(colData(sce_sim)$Cell), cluster=sce_sim$cluster)
}

clusters_annotated <- data.frame(sce@metadata$param_df,
                                 dimreduce_method=dimreduce_method,
                                 clustering_method=clustering_method,
                                 frac_genes=frac_genes,
                                 max_genes=max_genes,
                                 fc_percentile=fc_percentile,
                                 expr_percentile=expr_percentile,
                                 num_markers=n_markers,
                                 test_proportion=test_proportion,
                                 marker_setting=marker_setting,
                                 cluster_df,
                                 stringsAsFactors = FALSE)


## Write outputs

write.table(clusters_annotated, file = args$outfname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



