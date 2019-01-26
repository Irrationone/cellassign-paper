# Follicular evaluation plots for CellAssign

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

parser <- ArgumentParser(description = "Create follicular evaluation plot for CellAssign.")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Annotated SCE of follicular (subsetted for T cells, ideally)")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce <- readRDS(args$sce)
categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

nonother_types <- sort(setdiff(unique(sce$celltype), "other"))
sce <- sce %>%
  scater::mutate(celltype=factor(plyr::mapvalues(celltype, from = c("other"),
                                                 to = c("Unassigned")),
                                 levels = c(nonother_types, "Unassigned")))

# SC3 (markers) not shown because only one cluster is predicted
annotation_labels <- c("celltype", 
                       "all_sc3_cluster",
                       "all_seurat_0.8_cluster",
                       "all_subset_seurat_0.8_cluster",
                       "all_seurat_1.2_cluster",
                       "all_subset_seurat_1.2_cluster")
plot_titles <- c("CellAssign",
                 "SC3 (full)",
                 "Seurat (res = 0.8, full)",
                 "Seurat (res = 0.8, markers)",
                 "Seurat (res = 1.2, full)",
                 "Seurat (res = 1.2, markers)")
legend_titles <- c("CellAssign",
                   "cluster",
                   "cluster",
                   "cluster",
                   "cluster",
                   "cluster")

seurat_cluster_levels <- sort(unique(unlist(lapply(annotation_labels[3:length(annotation_labels)], function(x) as.character(unique(colData(sce)[,x]))))))
seurat_cluster_palette <- iwanthue(length(seurat_cluster_levels))
names(seurat_cluster_palette) <- seurat_cluster_levels

sc3_cluster_levels <- as.character(sort(as.numeric(unique(unlist(lapply(annotation_labels[2], function(x) as.character(unique(colData(sce)[,x]))))))))
sc3_cluster_palette <- iwanthue(length(sc3_cluster_levels))
names(sc3_cluster_palette) <- sc3_cluster_levels

palettes <- c(rep(list(categorical_palettes$celltype),
                  1),
              list(sc3_cluster_palette),
              rep(list(seurat_cluster_palette),
                  4))

follicular_celltype_plots <- lapply(seq_along(annotation_labels), function(i) {
  lab_col <- annotation_labels[i]
  colData(sce)[,lab_col] <- factor(colData(sce)[,lab_col])
  
  p <- plotReducedDim(sce,
                      use_dimred = args$dimreduce_type,
                      colour_by = lab_col,
                      point_alpha = 0.4, 
                      add_ticks = FALSE,
                      point_size = 2)
  p <- p + 
    guides(colour = FALSE,
           shape = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() + 
    guides(fill = guide_legend(title = legend_titles[i], 
                               nrow = max(2, floor(length(unique(sce$lab_col))/7)),
                               override.aes = list(alpha = 1, size = 2))) + 
    ggtitle(plot_titles[i])
  
  if (!is.na(palettes[[i]])) {
    p <- p + scale_fill_manual(values = palettes[[i]])
  }
  
  return(p)
})

## Marker plots

# Plots of marker gene expression
marker_genes <- c("CD8A", "ICOS", "CXCR5")
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
                      use_dimred = args$dimreduce_type,
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.4,
                      point_size = 1.5,
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

## Legend
marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Log normalized counts",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)

# Final plot

marker_plot_group <- cowplot::plot_grid(plotlist = marker_plots,
                                        ncol = 3,
                                        nrow = 1)

follicular_celltype_group <- cowplot::plot_grid(plotlist = follicular_celltype_plots,
                                          labels = letters[(1:length(follicular_celltype_plots))+1],
                                          nrow = 2,
                                          ncol = 3)

final_plot <- cowplot::plot_grid(marker_plot_group,
                                 marker_legend,
                                 follicular_celltype_group,
                                 labels = c('a', '', ''),
                                 nrow = 3,
                                 ncol = 1,
                                 rel_heights = c(0.5, 0.1, 1))


# Plot final plot
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


