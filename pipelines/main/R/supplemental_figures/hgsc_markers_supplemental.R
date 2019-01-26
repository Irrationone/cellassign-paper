# Supplemental marker figure for HGSC

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)
library(scales)
library(ggrastr)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create supplemental marker figure for HGSC")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

categorical_palettes <- cat_palettes()
gradient_colours <- scrna_expression_gradient()

# Marker plots

marker_genes <- c("PECAM1", "CD3D", "CD79A", 
                  "KLHDC8A", "ACTA2", "MYH11", 
                  "MCAM", "CD24", "WT1")

exprs <- logcounts(sce)[cellassign.utils::get_ensembl_id(marker_genes, sce),]
expr_limits <- c(min(exprs), max(exprs))

sce_tmp <- sce

if (!is.null(args$winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$winsorized_expression_threshold)
  expr_limits[2] <- min(expr_limits[2], args$winsorized_expression_threshold)
}

marker_plots <- lapply(marker_genes, function(mgene) {
  p <- plotReducedDim(sce_tmp,
                      use_dimred = args$dimreduce_type,
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.4,
                      point_size = 0.75,
                      add_ticks = FALSE)
  p$layers[[1]]$aes_params$colour <- NULL
  p$layers[[1]]$aes_params$shape <- 16
  p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$fill
  
  p <- p + 
    guides(fill = FALSE,
           colour = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})
names(marker_plots) <- marker_genes


## Expression values

marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Log normalized counts",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)

# Assemble plot

celltype_marker_plots_combined <- cowplot::plot_grid(
  plotlist = marker_plots[1:7],
  ncol = 4
)

cancer_marker_plots_combined <- cowplot::plot_grid(
  plotlist = marker_plots[8:9],
  ncol = 4
)

final_plot <- cowplot::plot_grid(celltype_marker_plots_combined,
                                 cancer_marker_plots_combined,
                                 marker_legend,
                                 labels = c('a', 'b', ''),
                                 nrow = 3,
                                 rel_heights = c(1, 0.5, 0.1))



# Plot to output file
pdf(args$outfname, width = 10, height = 8, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")








