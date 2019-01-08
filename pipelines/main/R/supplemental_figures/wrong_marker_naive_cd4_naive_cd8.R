# Supplemental wrong marker figure

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

parser <- ArgumentParser(description = "Create supplemental wrong marker figure.")
parser$add_argument('--wrongmarker_result_dir', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker")
parser$add_argument('--wrongmarker_result_dir2', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker (with 20 genes per type)")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

wrongmarker_result_dir <- args$wrongmarker_result_dir
wrongmarker_result_dir2 <- args$wrongmarker_result_dir2

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$clustering_methods

# Wrong marker figure (1)
wm_eval_measures <- load_annotation_files(wrongmarker_result_dir, pattern = "*_eval_measures.tsv")
wm_delta_vals <- load_annotation_files(wrongmarker_result_dir, pattern = "*_delta_compare.tsv")

wm_plots <- plot_simulation_performance(wm_eval_measures,
                                        measures = c("micro_f1",
                                                     "accuracy"),
                                        display_measure_names = c("F1",
                                                                  "Accuracy"),
                                        x_var = "wrong_marker_proportion")

wm_plot_cellassign <- wm_plots$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = guide_legend(title = "Method")) + 
  scale_fill_manual(values = clust_methods_palette)

# Wrong marker figure (2)
wm_eval_measures2 <- load_annotation_files(wrongmarker_result_dir2, pattern = "*_eval_measures.tsv")
wm_delta_vals2 <- load_annotation_files(wrongmarker_result_dir2, pattern = "*_delta_compare.tsv")

wm_plots2 <- plot_simulation_performance(wm_eval_measures2,
                                         measures = c("micro_f1",
                                                      "accuracy"),
                                         display_measure_names = c("F1",
                                                                   "Accuracy"),
                                         x_var = "wrong_marker_proportion")

wm_plot_cellassign2 <- wm_plots2$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = guide_legend(title = "Method")) + 
  scale_fill_manual(values = clust_methods_palette)

# Plot final plot
final_plot <- cowplot::plot_grid(wm_plot_cellassign,
                                 wm_plot_cellassign2,
                                 nrow = 2,
                                 rel_heights = c(0.5, 0.5), 
                                 labels = c('a', 'b'))

pdf(args$outfname, width = 7, height = 8, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


