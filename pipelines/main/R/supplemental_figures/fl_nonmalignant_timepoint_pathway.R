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
pval_colours <- pval_colour_gradient()

nonmalignant_classes <- c("cytotoxic", "helper", "follicular_helper", "b")
patients <- unique(tools::file_path_sans_ext(basename(de_timepoint_files)))

de_result <- lapply(patients, function(patient) {
  de_files <- de_timepoint_files[str_detect(de_timepoint_files, patient)]
  
  classes <- basename(dirname(de_files))
  classes_exist <- intersect(nonmalignant_classes, classes)
  
  res <- lapply(classes_exist, function(class) {
    readRDS(de_files[classes == class][1])
  })
  names(res) <- classes_exist
  return(res)
})
names(de_result) <- patients

pathway_tables <- list(
  de_result$FL1018$cytotoxic$pathway$up,
  de_result$FL2001$cytotoxic$pathway$up
)

legend_vals <- lapply(pathway_tables, function(x) {
  res <- x@result %>%
    dplyr::filter(p.adjust <= 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    head(30)
  return(list(padj=res$p.adjust, size=res$Count))
})

padj_limits <- c(min(sapply(legend_vals, function(x) x$padj)),
                 max(sapply(legend_vals, function(x) x$padj)))
size_limits <- c(min(sapply(legend_vals, function(x) x$size)),
                 max(sapply(legend_vals, function(x) x$size)))

fl1018_cytotoxic_up_plot <- plot_pathway_network(de_result$FL1018$cytotoxic$pathway$up %>%
                                                   filter_pathway_result(), shorten_names = TRUE,
                                                 guides = FALSE, colour_limits = padj_limits,
                                                 size_limits = size_limits, pval_colours = rev(pval_colours)) + 
  ggtitle("FL1018, CD8 T cells")
# fl1018_tfh_up_plot <- plot_pathway_network(de_result$FL1018$follicular_helper$pathway$up %>%
#                                              filter_pathway_result(), shorten_names = TRUE) + 
#   ggtitle("FL1018, Tfh")
# fl1018_cd4_up_plot <- plot_pathway_network(de_result$FL1018$helper$pathway$up %>%
#                                              filter_pathway_result(), shorten_names = TRUE) + 
#   ggtitle("FL1018, CD4 T cells")
# No pathways
# fl1018_b_up_plot <- plot_pathway_network(de_result$FL1018$b$pathway$up %>%
#                                            filter_pathway_result(), shorten_names = TRUE)

fl2001_cytotoxic_up_plot <- plot_pathway_network(de_result$FL2001$cytotoxic$pathway$up %>%
                                                   filter_pathway_result(), shorten_names = TRUE,
                                                 guides = FALSE, colour_limits = padj_limits,
                                                 size_limits = size_limits, pval_colours = rev(pval_colours)) + 
  ggtitle("FL2001, CD8 T cells") 
# fl2001_tfh_up_plot <- plot_pathway_network(de_result$FL2001$follicular_helper$pathway$up %>%
#                                              filter_pathway_result(), shorten_names = TRUE,
#                                            guides = FALSE) + 
#   ggtitle("FL2001, Tfh")
# # No pathways
# # fl2001_cd4_up_plot <- plot_pathway_network(de_result$FL2001$helper$pathway$up %>%
# #                                              filter_pathway_result(), shorten_names = TRUE)
# fl2001_b_up_plot <- plot_pathway_network(de_result$FL2001$b$pathway$up %>%
#                                            filter_pathway_result(), shorten_names = TRUE,
#                                          guides = FALSE) + 
#   ggtitle("FL2001, nonmalignant B")

# Legends

network_pval_legend <- cellassign.utils::ggsimplelegend(padj_limits,
                                                        colour_mapping = rev(pval_colours),
                                                        legend_title = "Adjusted p value",
                                                        type = "continuous") + 
  theme(legend.key.width = unit(2, "lines")) + 
  scale_fill_gradientn(colours = rev(pval_colours), limits = padj_limits, trans = "log10")
network_pval_legend <- cellassign.utils::extract_legend(network_pval_legend)

network_size_legend <- plot_pathway_network(de_result$FL1018$cytotoxic$pathway$up %>%
                                              filter_pathway_result(), shorten_names = TRUE,
                                            guides = FALSE, colour_limits = padj_limits,
                                            size_limits = size_limits, pval_colours = rev(pval_colours)) + 
  guides(colour = FALSE,
         size = guide_legend(title = "Genes", 
                             title.position = "top", 
                             title.hjust = 0.5))  + 
  scale_colour_gradientn(colours = rev(pval_colours),
                         trans = "log10",
                         limits = padj_limits)  + 
  scale_size_continuous(limits = size_limits, 
                        breaks = c(2, 5, 10)) + 
  theme_nature(fontsize = 7) + 
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.background = element_rect(fill = "transparent", 
                                         colour = NA), 
        legend.justification = "center") 
network_size_legend <- cellassign.utils::extract_legend(network_size_legend)


# Build combined plot

network_legend <- cowplot::plot_grid(network_pval_legend,
                                     network_size_legend,
                                     ncol = 2,
                                     rel_widths = c(0.5, 0.5))


final_plot <- cowplot::plot_grid(fl1018_cytotoxic_up_plot,
                                 network_legend,
                                 fl2001_cytotoxic_up_plot,
                                 labels = c('a', '', 'b'), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(1, .1, 1))

# Plot final plot
pdf(args$outfname, width = 7, height = 10)
plot(final_plot)
dev.off()

cat("Completed.\n")


