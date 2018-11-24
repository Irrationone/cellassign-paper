# Supplemental figure looking at kappa and lambda light chain expression

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

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create IGKC/IGLC figure for FL")
parser$add_argument('--sce_nonmalignant_bcell', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--cellassign_lambda_kappa', metavar='FILE', type='character',
                    help="Path to CellAssign lambda/kappa results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce_nonmalignant_bcell
sce <- readRDS(sce_path)
cellassign_lambda_kappa_results <- readRDS(args$cellassign_lambda_kappa)

categorical_palettes <- cat_palettes()

sce$class <- cellassign_lambda_kappa_results$cell_type

# IGKC/IGLC expression heatmap
marker_genes <- c("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7")

expression_heatmap <- plot_expression_heatmap(sce, 
                                              marker_genes = marker_genes,
                                              rowdat = rowData(sce),
                                              rowlabels = sce$class %>%
                                                plyr::mapvalues("other", "Unassigned"),
                                              n_sample = ncol(sce),
                                              label_name = "Class",
                                              annotation_colors = list(Class=categorical_palettes$light_chain_class))

# Plot final plot
pdf(expression_heatmap, width = 5, height = 5)
plot(expression_heatmap)
dev.off()

cat("Completed.\n")


