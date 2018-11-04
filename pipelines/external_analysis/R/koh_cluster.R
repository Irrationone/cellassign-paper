# Cellassign Koh data

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(Matrix)
library(Seurat)
library(SC3)
library(infotheo)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Cluster Koh data with cellassign.")
parser$add_argument('--sce', type = 'character', metavar = 'FILE',
                    help="Annotated SCE")
parser$add_argument('--rho', type = 'character', metavar='FILE',
                    help="Rho from CellAssign")
parser$add_argument('--clustering_methods', type = 'character', nargs ='+',
                    help="Clustering methods to run")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE")

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

sce_path <- args$sce
sce <- readRDS(sce_path)
rho <- read.table(args$rho, sep = "\t", row.names = 1, header = TRUE)
clustering_methods <- unlist(args$clustering_methods)

markers_to_use <- rownames(rho)
print(markers_to_use)

cluster_res <- lapply(clustering_methods, function(method) {
  sce <- sce %>%
    scater::mutate(Group = celltype,
                   Batch = 1)
  rowData(sce)$Gene <- rownames(rowData(sce))
  rowData(sce)$symbol <- rownames(rowData(sce))
  
  sce_full <- cluster_wrapper(sce, gene_subset = NULL, dimreduce_method = 'PCA', clustering_method = method, 
                              conda_env = 'r-tensorflow')
  
  sce_markers <- cluster_wrapper(sce, gene_subset = markers_to_use, dimreduce_method = 'PCA', 
                                 clustering_method = method, 
                                 conda_env = 'r-tensorflow')
  
  sce_full <- map_clusters(sce_full, method = "correlation", min_correlation = 0)
  corr_evaluation_measures_full <- compute_evaluation_measures(sce_full, 
                                                               truth_labels = sce_full$Group,
                                                               inferred_labels = sce_full$inferred_group)
  
  sce_markers <- map_clusters(sce_markers, method = "correlation", min_correlation = 0)
  corr_evaluation_measures_markers <- compute_evaluation_measures(sce_markers, 
                                                                  truth_labels = sce_markers$Group,
                                                                  inferred_labels = sce_markers$inferred_group)
  
  eval_measures <- rbind(data.frame(gene_set='full', corr_evaluation_measures_full),
                         data.frame(gene_set = 'markers', corr_evaluation_measures_markers))
  
  eval_measures <- data.frame(clustering_method=method, eval_measures)
  
  res <- list(clusters_full=sce_full$cluster, 
              clusters_markers=sce_markers$cluster,
              eval_measures=eval_measures)
  return(res)
})
names(cluster_res) <- clustering_methods

eval_measures_combined <- plyr::rbind.fill(lapply(cluster_res, function(x) x$eval_measures))

cluster_labels <- do.call('cbind', lapply(seq_along(cluster_res), function(i) {
  method <- names(cluster_res)[i]
  x <- cluster_res[[i]]
  df <- data.frame(x$clusters_full,
                   x$clusters_markers)
  colnames(df) <- paste0(method, c('_cluster_full', '_cluster_markers'))
  return(df)
}))


# Annotate
colData(sce) <- colData(sce) %>%
  data.frame(check.names = FALSE) %>%
  cbind(cluster_labels) %>%
  DataFrame(check.names = FALSE)


sce@metadata$evaluation_measures <- plyr::rbind.fill(sce@metadata$evaluation_measures,
                                                     eval_measures_combined)

# Save outputs
saveRDS(sce, file = args$outfname)

cat("Completed.\n")

