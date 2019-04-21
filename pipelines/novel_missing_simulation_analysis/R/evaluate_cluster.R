
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
parser$add_argument('--method', type='character',
                    help="Clustering method")
parser$add_argument('--marker_types', type='character', nargs = '+',
                    help="Marker cell types")
parser$add_argument('--outfname', metavar='FILE', type = 'character',
                    help="Output measures for clustering.")
parser$add_argument('--conda_env', type = 'character',
                    help="Conda environment", default = 'r-tensorflow')
args <- parser$parse_args()

marker_types <- unlist(args$marker_types)

sce <- readRDS(args$sce)
cluster_info <- fread(args$cluster_file)

clusters <- cluster_info %>%
  dplyr::select(cell, cluster) %>%
  dplyr::rename(Cell=cell)

sce_filtered <- sce %>%
  scater::filter(Cell %in% clusters$Cell)
colData(sce_filtered) <- colData(sce_filtered) %>%
  as.data.frame %>%
  dplyr::left_join(clusters) %>%
  DataFrame()

sce_filtered$Group <- as.character(sce_filtered$Group)
sce_filtered$Group[!sce_filtered$Group %in% marker_types] <- "unknown"

evaluation_measures <- compute_evaluation_measures(sce_filtered, 
                                                   truth_labels = sce_filtered$Group,
                                                   inferred_labels = sce_filtered$cluster)


run_params <- as.data.frame(cluster_info)[,which(!colnames(cluster_info) %in% c("cell", "cluster"))] %>%
  unique

eval_df <- data.frame(run_params,
                      evaluation_measures)

# Write outputs
write.table(eval_df, args$outfname, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




