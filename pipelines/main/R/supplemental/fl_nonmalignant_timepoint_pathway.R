# Plots of genes differentially expressed between timepoints for each nonmalignant cell type
# Pathway plots (ReactomePA)

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(Matrix)
library(org.Hs.eg.db)
library(ReactomePA)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create pathway plots for nonmalignant timepoint DE")
parser$add_argument('--de_timepoint_dir', type='character',
                    help="DE results for timepoint")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

de_timepoint_dir <- args$de_timepoint_dir
de_timepoint_files <- Sys.glob(file.path(de_timepoint_dir, "*", "*"))

categorical_palettes <- cat_palettes()




# Build combined plot

final_plot <- cowplot::plot_grid(plotlist = de_plots,
                                 labels = c('a', 'b'), 
                                 ncol = 1, 
                                 nrow = 2,
                                 rel_heights = c(0.5, 0.5))

# Plot final plot
pdf(args$outfname, width = 7, height = 10)
plot(final_plot)
dev.off()

cat("Completed.\n")


