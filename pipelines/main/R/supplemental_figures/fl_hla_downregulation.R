# Plot of HLA downregulation 

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
library(ggrepel)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create gene DE plots for HLA")
parser$add_argument('--de_timepoint_dir', type='character',
                    help="DE results for timepoint")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

de_timepoint_dir <- args$de_timepoint_dir
de_timepoint_files_malignant <- Sys.glob(file.path(de_timepoint_dir, "malignant", "*"))

highlight_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DRA", "HLA-DRB1", "CD74")
categorical_palettes <- cat_palettes()

patients <- tools::file_path_sans_ext(basename(de_timepoint_files_malignant))

de_plots <- lapply(seq_along(patients), function(i) {
  de_result <- readRDS(de_timepoint_files_malignant[i])
  
  gene_table <- de_result$gene %>%
    dplyr::mutate(Significance=plyr::mapvalues(is_significant,
                                               c(FALSE, TRUE), 
                                               c("P > 0.05", "P <= 0.05")))
  
  xlims <- c(-1,1) * max(abs(gene_table$logFC))
  
  p <- plot_markers_v2(gene_table,
                       top_n = 0,
                       top_direction = "equal",
                       highlight_genes = highlight_genes,
                       colour_col = "Significance", 
                       subtitle = "",
                       label_size = 2.5, 
                       point_alpha = 0.3,
                       force = 10,
                       min.segment.length = unit(0.1, 'lines')) + 
    scale_color_manual(values = categorical_palettes$significance) +
    guides(colour = guide_legend(title = "Significance",
                                 override.aes = list(alpha = 1))) + 
    ggtitle(patients[i]) + 
    scale_x_continuous(limits = xlims)
  
  return(p)
})
names(de_plots) <- patients


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


