# Overview figure for simulation results

# Figure comparing T1 to T2 in FL samples
# Components: cell cycle, differential expression, cell composition

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

parser <- ArgumentParser(description = "Create simulation figure.")
parser$add_argument('--deprob_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for DE prob")
parser$add_argument('--wrongmarker_result_dir', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir
wrongmarker_result_dir <- args$wrongmarker_result_dir

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

# DE prob figure
de_eval_measures <- load_annotation_files(deprob_result_dir, pattern = "*_eval_measures.tsv")
de_deltas <- load_annotation_files(deprob_result_dir, pattern = "*_delta_compare.tsv")

## TODO: When reruns have been done with ARI and NMI, add those too
de_plots <- plot_simulation_performance(de_eval_measures %>%
                                          dplyr::mutate(clustering_method=factor(clustering_method, 
                                                                                 levels =factor_orderings$clustering_methods)), 
                                        measures = c("v_measure",
                                                     "accuracy"),
                                        display_measure_names = c("V-measure",
                                                                  "Accuracy"),
                                        x_var = "de_prob")

## TODO: Check SC3 and gaussian -- code may not be working properly for it
de_plot_markers <- de_plots$markers + 
  guides(fill = FALSE) + 
  xlab("% DE per group") + 
  scale_fill_manual(values = categorical_palettes$clustering_methods)
  
de_plot_full <- de_plots$full + 
  guides(fill = FALSE) + 
  xlab("% DE per group") + 
  scale_fill_manual(values = categorical_palettes$clustering_methods)


de_plot_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$clustering_methods),
                                                   colour_mapping = unname(categorical_palettes$clustering_methods),
                                                   legend_title = "Method", legend_rows = 2, fontsize = 7)
de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)

# Wrong marker figure
wm_eval_measures <- load_annotation_files(wrongmarker_result_dir, pattern = "*_eval_measures.tsv")
wm_delta_vals <- load_annotation_files(wrongmarker_result_dir, pattern = "*_delta_compare.tsv")

wm_plots <- plot_simulation_performance(wm_eval_measures, 
                                        measures = c("v_measure",
                                                     "accuracy"),
                                        display_measure_names = c("V-measure",
                                                                  "Accuracy"),
                                        x_var = "wrong_marker_proportion")

wm_plot_cellassign <- wm_plots$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = guide_legend(title = "Method"))


# Final plot
de_plots_labeled <- cowplot::plot_grid(de_plot_markers, de_plot_full, de_plot_legend,
                                       labels = c('a', 'b', ''),
                                       ncol = 1,
                                       nrow = 3,
                                       rel_heights = c(1, 1, 0.2))

final_plot <- cowplot::plot_grid(de_plots_labeled, wm_plot_cellassign, 
                                 labels = c('', 'c'), 
                                 ncol = 1, 
                                 nrow = 2,
                                 rel_heights = c(0.67, 0.33))


# Plot final plot
pdf(args$outfname, width = 10, height = 10)
plot(final_plot)
dev.off()

cat("Completed.\n")


