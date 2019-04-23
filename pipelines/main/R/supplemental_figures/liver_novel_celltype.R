# Liver novel celltype plots

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

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Liver novel celltype plots.")
parser$add_argument('--sce_liver', metavar='FILE', type='character',
                    help="SCE of liver")
parser$add_argument('--marker_types', type='character', nargs = '+',
                    help="Celltypes with marker genes")
parser$add_argument('--celltypes_novel1', type='character', nargs = '+',
                    help="Celltypes in data (novel 1)")
parser$add_argument('--celltypes_novel8', type='character', nargs = '+',
                    help="Celltypes in data (novel 8)")
parser$add_argument('--cellassign_fit_novel1', metavar='FILE', type='character',
                    help="CellAssign fit to 3 cell types")
parser$add_argument('--cellassign_fit_novel8_raw', metavar='FILE', type='character',
                    help="CellAssign fit to 4 cell types")
parser$add_argument('--cellassign_fit_novel8', metavar='FILE', type='character',
                    help="CellAssign fit to 4 cell types")
parser$add_argument('--scina_fit_novel1', metavar='FILE', type='character',
                    help="SCINA fit to 3 cell types")
parser$add_argument('--scina_fit_novel8_raw', metavar='FILE', type='character',
                    help="SCINA fit to 4 cell types")
parser$add_argument('--scina_fit_novel8_revised', metavar='FILE', type='character',
                    help="SCINA fit to 11 cell types")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_liver <- readRDS(args$sce_liver)
marker_types <- unlist(args$marker_types)

model_fits <- list(
  'cellassign_fit_novel1'=readRDS(args$cellassign_fit_novel1),
  'cellassign_fit_novel8_raw'=readRDS(args$cellassign_fit_novel8_raw),
  'cellassign_fit_novel8'=readRDS(args$cellassign_fit_novel8),
  'scina_fit_novel1'=readRDS(args$scina_fit_novel1),
  'scina_fit_novel8_raw'=readRDS(args$scina_fit_novel8_raw),
  'scina_fit_novel8'=readRDS(args$scina_fit_novel8)
)

celltypes_used <- list(
  'cellassign_fit_novel1'=unlist(args$celltypes_novel1),
  'cellassign_fit_novel8_raw'=unlist(args$celltypes_novel8),
  'cellassign_fit_novel8'=unlist(args$celltypes_novel8),
  'scina_fit_novel1'=unlist(args$celltypes_novel1),
  'scina_fit_novel8_raw'=unlist(args$celltypes_novel8),
  'scina_fit_novel8'=unlist(args$celltypes_novel8)
)

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

liver_celltype_palette <- categorical_palettes$liver_celltypes
liver_celltype_palette <- c(liver_celltype_palette, 'Unassigned'='gray60')

plot_titles <- c("CellAssign (+ NK cells)",
                 "CellAssign (all cells)",
                 "CellAssign (all cells, revised markers)",
                 "SCINA (+ NK cells)",
                 "SCINA (all cells)",
                 "SCINA (all cells, revised markers)")

celltype_plots <- lapply(seq_along(names(model_fits)), function(i) {
  fit_name <- names(model_fits)[i]
  fit <- model_fits[[fit_name]]
  celltype_subset <- celltypes_used[[fit_name]]
  
  if ("cell_type" %in% names(fit)) {
    celltypes <- fit$cell_type %>%
      plyr::mapvalues("other", "Unassigned")
  } else if ("cell_labels" %in% names(fit)) {
    celltypes <- fit$cell_labels %>%
      plyr::mapvalues("unknown", "Unassigned")
  } else {
    stop("Cannot find cell labels.")
  }
  
  sce_liver_subset <- sce_liver %>%
    scater::filter(celltype %in% celltype_subset)
  
  celltype_levels <- sort(unique(c(sce_liver_subset$celltype, celltypes)))
  
  cont_table <- table(factor(sce_liver_subset$celltype, levels = celltype_levels), 
                      factor(celltypes, levels = celltype_levels))
  acc <- sum(diag(cont_table))/sum(cont_table)
  micro_f1 <- microF1(cont_table)
  
  sce_liver_subset$assigned_cluster <- celltypes
  
  p <- plotReducedDim(sce_liver_subset,
                      use_dimred = args$dimreduce_type,
                      colour_by = "assigned_cluster",
                      point_alpha = 0.4, 
                      point_size = 1.5)
  p <- p + 
    guides(colour = FALSE,
           shape = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() + 
    scale_fill_manual(values = liver_celltype_palette) + 
    guides(fill = FALSE) + 
    ggtitle(plot_titles[i])
  
  plot_label <- paste0("Accuracy = ",  format(acc, digits = 3), "\n", "F1 = ", format(micro_f1, digits = 3))
  p <- p + annotate(geom = 'text', x = Inf, y = Inf, hjust = 1, vjust = 1.5, label = plot_label, parse = FALSE,
                    size = 2.5)
  
  return(p)
})
names(celltype_plots) <- names(model_fits)

master_label_plot <- plotReducedDim(sce_liver,
                                    use_dimred = args$dimreduce_type,
                                    colour_by = "celltype",
                                    point_alpha = 0.4, 
                                    point_size = 1.5) +
  guides(colour = FALSE,
         shape = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = liver_celltype_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Celltype")

celltype_legend <- cellassign.utils::ggsimplelegend(names(liver_celltype_palette),
                                                    colour_mapping = unname(liver_celltype_palette),
                                                    legend_title = "Celltype", legend_rows = 3, fontsize = 7)
celltype_legend <- cellassign.utils::extract_legend(celltype_legend)

# Final plot
complete_plotlist <- c(list(master_label_plot), celltype_plots)
combined_plots <- cowplot::plot_grid(plotlist = complete_plotlist,
                                     labels = letters[1:length(complete_plotlist)],
                                     nrow = 3,
                                     ncol = 3)

final_plot <- cowplot::plot_grid(combined_plots,
                                 celltype_legend,
                                 nrow = 2,
                                 rel_heights = c(0.9, 0.1))



# Plot final plot
pdf(args$outfname, width = 10, height = 14, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


