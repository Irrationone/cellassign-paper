
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

library(scrna.utils)
library(cellassign)
library(cellassign.utils)
library(scrna.sceutils)
library(Seurat)

library(methods)
library(argparse)

parser <- ArgumentParser(description = "Cluster SCE with CellAssign")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Input SCE object")
parser$add_argument('--cluster_file', metavar='FILE', type='character',
                    help="Cluster results")
parser$add_argument('--outfname', metavar='FILE', type = 'character',
                    help="Output measures for clustering.")
parser$add_argument('--conda_env', type = 'character',
                    help="Conda environment", default = 'r-tensorflow')
args <- parser$parse_args()


map_and_evaluate_clusters <- function(sce_sim, min_count = 0) {
  # Remove genes that are lowly expressed, and cells with no UMIs
  sce_sim <- sce_sim[rowSums(counts(sce_sim)) >= min_count,colSums(counts(sce_sim)) > 0]
  
  nclusts <- length(unique(sce_sim$cluster))
  
  if (nclusts > 1 & !all(is.na(sce_sim$cluster))) {
    sce_de <- map_clusters(sce_sim, method = "de2", FDR_cutoff = 0.05)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0.1)
  } else {
    max_group <- names(which.max(table(sce_sim$Group)))[1]
    
    sce_de <- sce %>% scater::mutate(inferred_group=max_group)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0.1)
  }
  
  de_evaluation_measures <- compute_evaluation_measures(sce_de, 
                                                        truth_labels = sce_de$Group,
                                                        inferred_labels = sce_de$inferred_group)
  corr_evaluation_measures <- compute_evaluation_measures(sce_correlation, 
                                                          truth_labels = sce_correlation$Group,
                                                          inferred_labels = sce_correlation$inferred_group)
  
  evaluation_measures <- rbind(data.frame(mapping_type = "de", de_evaluation_measures, stringsAsFactors = FALSE), data.frame(mapping_type = "correlation", corr_evaluation_measures, stringsAsFactors = FALSE))
  
  return(evaluation_measures)
}

sce <- readRDS(args$sce)
cluster_info <- fread(args$cluster_file)

clusters <- cluster_info %>%
  dplyr::select(Cell, cluster)

sce_filtered <- sce %>%
  scater::filter(Cell %in% clusters$Cell)
colData(sce_filtered) <- colData(sce_filtered) %>%
  as.data.frame %>%
  dplyr::left_join(clusters) %>%
  DataFrame()

if (all(sce_filtered$cluster %in% sce_filtered$Group)) {
  message("Evaluating directly ...")
  evaluation_measures <- compute_evaluation_measures(sce_filtered, 
                                                     truth_labels = sce_filtered$Group,
                                                     inferred_labels = sce_filtered$cluster)
} else {
  message("Mapping then evaluating ...")
  evaluation_measures <- map_and_evaluate_clusters(sce_filtered)
}

run_params <- cluster_info[,!colnames(cluster_info) %in% c("Cell", "cluster")] %>%
  unique

eval_df <- data.frame(run_params,
                      evaluation_measures)

# Write outputs
write.table(eval_df, args$outfname, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




