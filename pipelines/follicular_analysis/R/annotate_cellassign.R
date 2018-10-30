# Annotate FL SCE with CellAssign-based predictions
# NOTE: This script is FL-specific

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

parser <- ArgumentParser(description = "Annotate SCE with CellAssign results")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--broad', metavar='FILE', type='character',
                    help="Path to broad assignments")
parser$add_argument('--specific', metavar='FILE', type='character',
                    help="Path to specific assignments")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

broad_assignments <- readRDS(args$broad)
specific_assignments <- readRDS(args$specific)

# Add MAP assignments
sce$cellassign_cluster_broad <- broad_assignments$cell_type
sce$cellassign_cluster_specific <- specific_assignments$cell_type

# Add probabilities
broad_gamma <- as.data.frame(broad_assignments$mle_params$gamma)
colnames(broad_gamma) <- paste0(colnames(broad_gamma), " (broad)")
specific_gamma <- as.data.frame(specific_assignments$mle_params$gamma)

sce@colData <- bind_cols(sce@colData %>% as.data.frame, 
                         broad_gamma,
                         specific_gamma) %>% DataFrame(check.names = FALSE)

# Set master celltype label
sce$celltype <- sce$cellassign_cluster_specific

# Save annotated SCE
saveRDS(sce, args$outfname)

cat("Completed.\n")



