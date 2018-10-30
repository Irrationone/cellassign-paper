# Label FL cells as malignant or nonmalignant based on Ig-chain status
# In reality this needs to be updated EVERY TIME

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
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Extract B cells
b_cells <- sce %>% 
  scater::filter(celltype == "B cells")

# Cluster deterministically with Seurat
set.seed(19279)
b_cells <- cluster_wrapper(b_cells, 
                           gene_subset = NULL, 
                           dimreduce_method = "PCA", 
                           clustering_method = "seurat",
                           seurat_resolution = 0.8)

nonmalignant_clusters <- c("6")

B_mapping <- colData(b_cells) %>%
  data.frame(check.names = FALSE) %>%
  dplyr::select(sample_barcode, cluster) %>%
  dplyr::mutate(malignant_status_manual=ifelse(cluster %in% nonmalignant_clusters, "nonmalignant", "malignant")) %>%
  dplyr::select(-c(cluster))

colData(sce) <- colData(sce) %>%
  data.frame(check.names = FALSE) %>%
  dplyr::left_join(B_mapping) %>%
  DataFrame

sce$malignant_status_manual[is.na(sce$malignant_status_manual)] <- "nonmalignant"

sce <- sce %>%
  scater::mutate(
    celltype_full=ifelse(malignant_status_manual == "malignant", "B cells (malignant)", as.character(celltype))
  )

# Save malignant status
saveRDS(sce, args$outfname)

cat("Completed.\n")



