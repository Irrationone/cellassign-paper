# HGSC evaluation plots for CellAssign

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

parser <- ArgumentParser(description = "Create HGSC evaluation plot for CellAssign.")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Annotated SCE of HGSC")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
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

sc3_special_cluster <- "40"

seurat_cluster_levels <- sort(unique(unlist(lapply(annotation_labels[3:length(annotation_labels)], function(x) as.character(unique(colData(sce)[,x]))))))
seurat_cluster_palette <- iwanthue(length(seurat_cluster_levels))
names(seurat_cluster_palette) <- seurat_cluster_levels

sc3_cluster_levels <- as.character(sort(as.numeric(unique(unlist(lapply(annotation_labels[2], function(x) as.character(unique(colData(sce)[,x]))))))))
sc3_cluster_palette <- c(viridisLite::viridis(n = length(sc3_cluster_levels)-1), "#D54E37")
names(sc3_cluster_palette) <- sc3_cluster_levels

palettes <- c(rep(list(categorical_palettes$hgsc_celltype),
                  1),
              list(sc3_cluster_palette),
              rep(list(seurat_cluster_palette),
                  4))

hgsc_celltype_plots <- lapply(seq_along(annotation_labels), function(i) {
  lab_col <- annotation_labels[i]
  
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
                               override.aes = list(alpha = 1))) + 
    ggtitle(plot_titles[i])
  
  if (!is.na(palettes[[i]])) {
    p <- p + scale_fill_manual(values = palettes[[i]])
  }
  
  return(p)
})




# Final plot
final_plot <- cowplot::plot_grid(plotlist = hgsc_celltype_plots,
                                 labels = letters[1:length(hgsc_celltype_plots)],
                                 nrow = 3,
                                 ncol = 2,
                                 rel_widths = rep(1,2),
                                 rel_heights = rep(1,3))


# Plot final plot
pdf(args$outfname, width = 7, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


