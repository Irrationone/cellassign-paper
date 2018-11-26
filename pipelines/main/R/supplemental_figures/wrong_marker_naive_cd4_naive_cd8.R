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
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

wrongmarker_result_dir <- args$wrongmarker_result_dir

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$clustering_methods

# Wrong marker figure
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

# Plot final plot
pdf(args$outfname, width = 7, height = 5, useDingbats = FALSE)
plot(wm_plot_cellassign)
dev.off()

cat("Completed.\n")


