# Liver all celltype plots

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

parser <- ArgumentParser(description = "Liver all celltype plots.")
parser$add_argument('--sce_liver', metavar='FILE', type='character',
                    help="SCE of liver")
parser$add_argument('--cellassign_fit_raw', metavar='FILE', type='character',
                    help="CellAssign fit with raw markers")
parser$add_argument('--cellassign_fit_revised', metavar='FILE', type='character',
                    help="CellAssign fit with revised markers")
parser$add_argument('--scina_fit_raw', metavar='FILE', type='character',
                    help="SCINA fit with raw markers")
parser$add_argument('--scina_fit_revised', metavar='FILE', type='character',
                    help="SCINA fit with revised markers")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_liver <- readRDS(args$sce_liver)

model_fits <- list(
  'cellassign_fit_raw'=readRDS(args$cellassign_fit_raw),
  'cellassign_fit_revised'=readRDS(args$cellassign_fit_revised),
  'scina_fit_raw'=readRDS(args$scina_fit_raw),
  'scina_fit_revised'=readRDS(args$scina_fit_revised)
)

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

liver_celltype_palette <- categorical_palettes$liver_celltypes
liver_celltype_palette <- c(liver_celltype_palette, 'Unassigned'='gray60')

plot_titles <- c("CellAssign",
                 "CellAssign (revised markers)",
                 "SCINA",
                 "SCINA (revised markers)")

celltype_plots <- lapply(seq_along(names(model_fits)), function(i) {
  fit_name <- names(model_fits)[i]
  fit <- model_fits[[fit_name]]
  
  if ("cell_type" %in% names(fit)) {
    celltypes <- fit$cell_type %>%
      plyr::mapvalues("other", "Unassigned")
  } else if ("cell_labels" %in% names(fit)) {
    celltypes <- fit$cell_labels %>%
      plyr::mapvalues("unknown", "Unassigned")
  } else {
    stop("Cannot find cell labels.")
  }
  
  celltype_levels <- sort(unique(c(sce_liver$celltype, celltypes)))
  
  cont_table <- table(factor(sce_liver$celltype, levels = celltype_levels), 
                      factor(celltypes, levels = celltype_levels))
  acc <- sum(diag(cont_table))/sum(cont_table)
  micro_f1 <- microF1(cont_table)
  
  sce_liver$assigned_cluster <- celltypes
  
  p <- plotReducedDim(sce_liver,
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
  
  plot_label <- paste0("Acc = ",  format(acc, digits = 3), "\n", "F1 = ", format(micro_f1, digits = 3))
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

combined_plots <- cowplot::plot_grid(plotlist = celltype_plots,
                                     labels = letters[2:(length(celltype_plots)+1)],
                                     nrow = 2,
                                     ncol = 2)

final_plot <- cowplot::plot_grid(cowplot::plot_grid(master_label_plot, labels = c('a', ''), ncol = 2),
                                 combined_plots,
                                 celltype_legend,
                                 nrow = 3,
                                 rel_heights = c(0.45, 0.9, 0.2))



# Plot final plot
pdf(args$outfname, width = 10, height = 12, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


