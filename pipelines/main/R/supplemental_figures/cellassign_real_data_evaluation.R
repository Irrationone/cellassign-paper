# Real data evaluation plots for CellAssign

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

parser <- ArgumentParser(description = "Create real data evaluation plot for CellAssign.")
parser$add_argument('--koh_annotated', metavar='FILE', type='character',
                    help="Annotated SCE of Koh")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

koh_annotated_path <- args$koh_annotated

sce_koh <- readRDS(koh_annotated_path)
categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

evaluation_measures <- sce_koh@metadata$evaluation_measures %>%
  dplyr::mutate(full_label=ifelse(is.na(gene_set),
                                  as.character(clustering_method),
                                  paste0(clustering_method, "_", gene_set)))

# SC3 (markers) not shown because only one cluster is predicted
annotation_labels <- c("celltype", 
                       "cellassign_cluster",
                       "SC3_cluster_full",
                       "seurat_0.8_cluster_full",
                       "seurat_0.8_cluster_markers",
                       "seurat_1.2_cluster_full",
                       "seurat_1.2_cluster_markers")
plot_titles <- c("Celltype",
                 "CellAssign",
                 "SC3 (full)",
                 "Seurat (res = 0.8, full)",
                 "Seurat (res = 0.8, markers)",
                 "Seurat (res = 1.2, full)",
                 "Seurat (res = 1.2, markers)")
legend_titles <- c("Celltype",
                   "CellAssign",
                   "cluster",
                   "cluster",
                   "cluster",
                   "cluster",
                   "cluster")
eval_measure_labels <- c(NA,
                         "cellassign-shrinkage",
                         "SC3_full",
                         "seurat_0.8_full",
                         "seurat_0.8_markers",
                         "seurat_1.2_full",
                         "seurat_1.2_markers")

cluster_levels <- sort(unique(unlist(lapply(annotation_labels[3:length(annotation_labels)], function(x) as.character(unique(colData(sce_koh)[,x]))))))
cluster_palette <- iwanthue(length(cluster_levels))
names(cluster_palette) <- cluster_levels

palettes <- c(rep(list(categorical_palettes$koh_celltype),
                  2),
              rep(list(cluster_palette),
                  5))

koh_celltype_plots <- lapply(seq_along(annotation_labels), function(i) {
  lab_col <- annotation_labels[i]
  eval_measure_label <- eval_measure_labels[i]
  
  p <- plotReducedDim(sce_koh,
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
    scale_fill_manual(values = palettes[[i]]) + 
    guides(fill = guide_legend(title = legend_titles[i], 
                               override.aes = list(alpha = 1))) + 
    ggtitle(plot_titles[i])
  
  if (!is.na(eval_measure_label)) {
    eval_scores <- evaluation_measures %>% 
      dplyr::filter(full_label == eval_measure_label)
    plot_label <- paste0("Accuracy = ",  format(eval_scores$accuracy, digits = 3), "\n", "F1 = ", format(eval_scores$macro_f1, digits = 3))
    
    p <- p + annotate(geom = 'text', x = Inf, y = Inf, hjust = 1, vjust = 1.5, label = plot_label, parse = FALSE,
                      size = 2.5)
  }
  
  return(p)
})




# Final plot
final_plot <- cowplot::plot_grid(plotlist = koh_celltype_plots,
                                 labels = letters[1:length(koh_celltype_plots)],
                                 nrow = 3,
                                 ncol = 3,
                                 rel_widths = rep(1,3),
                                 rel_heights = rep(1,3))


# Plot final plot
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


