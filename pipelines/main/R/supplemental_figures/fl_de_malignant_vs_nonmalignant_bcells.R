# Plot of genes upregulated in malignant vs. nonmalignant B cells (highlight BCL2, etc.)

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
library(ggrepel)
library(ReactomePA)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create malignant vs. nonmalignant supplemental figure")
parser$add_argument('--de_malignant_timepoint_dir', type='character',
                    help="DE results for malignant-timepoint interaction comparisons")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

de_malignant_timepoint_dir <- args$de_malignant_timepoint_dir
de_malignant_timepoint_files <- Sys.glob(file.path(de_malignant_timepoint_dir, "*", "*"))
highlight_genes <- c("BCL2", "BCL6")

categorical_palettes <- cat_palettes()

patients <- tools::file_path_sans_ext(basename(de_malignant_timepoint_files))

de_plots <- lapply(seq_along(patients), function(i) {
  de_result <- readRDS(de_malignant_timepoint_files[i])
  
  malignant_gene_table <- de_result$malignant$gene %>%
    dplyr::mutate(Significance=plyr::mapvalues(is_significant,
                                               c(FALSE, TRUE), 
                                               c("P > 0.05", "P <= 0.05")))
  
  xlims <- c(-1,1) * max(abs(malignant_gene_table$logFC))
  
  p <- plot_markers_v2(malignant_gene_table,
                       top_n = 5,
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


