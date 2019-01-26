# Figure showing expression of markers in RLN (using scvis space)

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)
library(ggrastr)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "RLN marker expression plot")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for RLN + FL (merged)")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

sce <- sce %>%
  scater::filter(set_name == "RLN")

# Plots of marker gene expression
marker_genes <- c("CD2", "CD3D", "CD3E", "CD79A", "MS4A1", "CD19")
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
                      use_dimred = "scvis",
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
    xlab("scvis-1") + 
    ylab("scvis-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})

# Legends

## Expression values

marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Log normalized counts",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)


# Final plot

marker_plot_group <- cowplot::plot_grid(plotlist = marker_plots,
                                        ncol = 3,
                                        nrow = 2)

final_plot <- cowplot::plot_grid(marker_plot_group,
                                 marker_legend,
                                 nrow = 2,
                                 ncol = 1,
                                 rel_heights = c(0.9, 0.1),
                                 labels = c('', ''))

pdf(args$outfname, width = 8, height = 4.5, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")

