# Overview figure showing FL cell types

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create overview figure for FL")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

categorical_palettes <- cat_palettes()

# Plot of progression
df <- data.frame(
  patient=c('FL1018', 'FL1018', 'FL1018', 'FL1018'),
  timepoint=c('P', 'T1', 'T2', 'D'),
  time=c(0, 3, 4.5, 9)
)

timepoint_plot <- ggplot(df, aes(x=time, y=patient)) + 
  geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) + 
  geom_point(aes(colour=timepoint), size = 5) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Time since diagnosis (years)") + 
  ylab("") + 
  scale_colour_manual(values = categorical_palettes$timepoint, limits = c('P', 'T1', 'T2')) + 
  guides(colour = guide_legend(title = "Timepoint")) + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = 'bold'))

# Plot of timepoint
dr_timepoint <- plotReducedDim(sce, use_dimred = "UMAP", colour_by = "timepoint", point_alpha = 0.5, add_ticks = FALSE)
dr_timepoint <- dr_timepoint + 
  geom_rug(alpha = 0.1, colour = "gray20") +
  guides(fill = guide_legend(title = "Timepoint")) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = categorical_palettes$timepoint)

# Plot of celltype assignments
nonother_types <- sort(setdiff(unique(sce$celltype_full), "other"))
dr_celltype <- plotReducedDim(sce %>%
                                scater::mutate(celltype_full=factor(plyr::mapvalues(celltype_full, from = c("other"),
                                                                                    to = c("Unassigned")),
                                                                    levels = c(nonother_types, "Unassigned"))), 
                              use_dimred = "UMAP",
                              colour_by = "celltype_full",
                              point_alpha = 0.5, 
                              add_ticks = FALSE)
dr_celltype <- dr_celltype + 
  geom_rug(alpha = 0.1, colour = "gray20") +
  guides(fill = guide_legend(title = "Predicted celltype")) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = categorical_palettes$celltype)

# Plots of marker gene expression
marker_genes <- c("CD79A", "CD3D", "IGKC", "IGLC2")
exprs <- logcounts(sce)[cellassign.utils::get_ensembl_id(marker_genes, sce),]
expr_limits <- c(min(exprs), max(exprs))

gradient_colours <- scrna_expression_gradient()

marker_plots <- lapply(marker_genes, function(mgene) {
  p <- plotReducedDim(sce,
                 use_dimred = "UMAP",
                 colour_by = cellassign.utils::get_ensembl_id(mgene, sce),
                 point_alpha = 0.5,
                 point_size = 0.5,
                 add_ticks = FALSE)
  p <- p + 
    guides(fill = FALSE) + 
    xlab("UMAP-1") + 
    ylab("UMAP-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_fill_gradientn(colours = gradient_colours, 
                         limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})

# Assemble plot

## DR plots
dr_plots <- cowplot::plot_grid(dr_timepoint, dr_celltype, ncol = 2, nrow = 1, labels = c('b', 'c'))
marker_gene_plots <- cowplot::plot_grid(plotlist = marker_plots, ncol = 4, nrow = 1, labels = c('d', 'e', 'f', 'g'))

final_plot <- cowplot::plot_grid(timepoint_plot, dr_plots, marker_gene_plots, 
                                 labels = c('a', '', ''), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(0.5, 1, 0.5))

# Plot to output file
pdf(args$outfname, width = 10, height = 10)
final_plot
dev.off()

cat("Completed.\n")



