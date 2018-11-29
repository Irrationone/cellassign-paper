# Benchmarking plot

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

parser <- ArgumentParser(description = "Benchmarking figure")
parser$add_argument('--benchmark_times', type = 'character', metavar='FILE',
                    help="Benchmarking timing results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

benchmark_file <- args$benchmark_times
timing_df <- fread(benchmark_file)

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$clustering_methods[deprob_methods]

## Create subfigure for num cells scaling

timing_df_ncells <- timing_df %>% 
  dplyr::filter(num_groups == 2,
                num_genes == 1e4)

ggplot(timing_df_ncells, aes(x=num_cells, y = time, fill=factor(max_genes))) + 
  geom_boxplot(aes(fill=factor(max_genes))) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Number of cells") + 
  ylab("Time to convergence (seconds)")

final_plot <- cowplot::plot_grid(de_plots_labeled, 
                                 bottom_row,
                                 delta_plot_legend,
                                 labels = c('', ''), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(0.67, 0.33, 0.05))


# Plot final plot
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


