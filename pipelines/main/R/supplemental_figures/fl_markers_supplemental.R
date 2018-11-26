# Figure showing expression of additional markers

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(grid)
library(ggplotify)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "FL supplemental marker expression plot")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--cellassign_results', type = 'character', metavar = 'FILE',
                    help="CellAssign specific assignments results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
specific_assignments <- readRDS(args$cellassign_results)

categorical_palettes <- cat_palettes()

# Plots of marker gene expression
marker_genes <- c("CD2", "MS4A1", "CD8A", "GZMA", "CD4", "CXCR5", "ICOS")
exprs <- logcounts(sce)[cellassign.utils::get_ensembl_id(marker_genes, sce),]
expr_limits <- c(min(exprs), max(exprs))

gradient_colours <- scrna_expression_gradient()

sce_tmp <- sce

if (!is.null(args$winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$winsorized_expression_threshold)
  expr_limits[2] <- min(expr_limits[2], args$winsorized_expression_threshold)
}

marker_plots <- lapply(marker_genes, function(mgene) {
  p <- plotReducedDim(sce_tmp,
                      use_dimred = "UMAP",
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.2,
                      point_size = 0.5,
                      add_ticks = FALSE)
  p$layers[[1]]$aes_params$colour <- NULL
  p$layers[[1]]$aes_params$shape <- 16
  p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$fill
  
  p <- p + 
    guides(fill = FALSE,
           colour = FALSE) + 
    xlab("UMAP-1") + 
    ylab("UMAP-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})

# Marker gene expression heatmap
celltype_labels <- sce$celltype_full %>%
  plyr::mapvalues(from = c("other"), to = c("Unassigned"))
expression_heatmap <- plot_expression_heatmap(sce, 
                                              marker_genes = rownames(specific_assignments$mle_params$delta),
                                              rowdat = rowData(sce),
                                              rowlabels = celltype_labels,
                                              n_sample = ncol(sce),
                                              label_name = "Celltype",
                                              annotation_colors = list(Celltype=categorical_palettes$celltype[unique(celltype_labels)]))

# Legends

## Expression values

marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Expression",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)


# Final plot

marker_plot_group <- cowplot::plot_grid(plotlist = marker_plots,
                                        ncol = 4,
                                        nrow = 2)

final_plot <- cowplot::plot_grid(marker_plot_group,
                                 marker_legend,
                                 as.grob(expression_heatmap),
                                 nrow = 3,
                                 ncol = 1,
                                 rel_heights = c(0.6, 0.05, 0.4),
                                 labels = c('a', '', 'b'))

pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")

