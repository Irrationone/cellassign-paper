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


## Create subfigure for num cells scaling

timing_df_ncells <- timing_df %>% 
  dplyr::filter(num_groups == 2,
                num_genes == 1e4)

cell_timing_plot <- ggplot(timing_df_ncells, aes(x=factor(num_cells), y = time, fill=factor(max_genes))) + 
  geom_boxplot(outlier.size = -1) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), alpha = 0.4, size = 1) +
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Number of cells") + 
  ylab("Time to convergence (seconds)") + 
  scale_fill_manual(values = categorical_palettes$num_markers) +
  guides(fill = FALSE)

timing_df_ngroups <- timing_df %>% 
  dplyr::filter(num_cells == 1000,
                num_genes == 1e4)

group_timing_plot <- ggplot(timing_df_ngroups, aes(x=factor(num_groups), y = time, fill=factor(max_genes))) + 
  geom_boxplot(width = 0.5, outlier.size = -1) + 
  geom_jitter(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1, jitter.height = 0), alpha = 0.4, size = 1) +
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Number of cell types") + 
  ylab("Time to convergence (seconds)") + 
  scale_fill_manual(values = categorical_palettes$num_markers) +
  guides(fill = FALSE)

## Legends

marker_legend <- cellassign.utils::ggsimplelegend(unique(timing_df_ngroups$max_genes),
                                                  colour_mapping = categorical_palettes$num_markers,
                                                  legend_title = "Markers per celltype",
                                                  type = "discrete") 
marker_legend <- cellassign.utils::extract_legend(marker_legend)

final_plot <- cowplot::plot_grid(cell_timing_plot, 
                                 marker_legend,
                                 group_timing_plot,
                                 labels = c('a', '', 'b'), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(1, 0.1, 1))


# Plot final plot
pdf(args$outfname, width = 8, height = 8, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


