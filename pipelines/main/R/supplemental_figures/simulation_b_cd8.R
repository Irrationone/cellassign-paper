# Simulation results for B vs. CD8

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
parser$add_argument('--delta_deprobs', type='double', nargs ='+',
                    help="DE probs to show delta plots for.")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir
delta_deprobs <- unlist(args$delta_deprobs)
deprob_methods <- unlist(args$deprob_methods)

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$clustering_methods[deprob_methods]

# DE prob figure
de_eval_measures <- load_annotation_files(deprob_result_dir, pattern = "*_eval_measures.tsv")
de_deltas <- load_annotation_files(deprob_result_dir, pattern = "*_delta_compare.tsv")

## Only use methods that were selected
de_eval_measures <- de_eval_measures %>%
  dplyr::filter(clustering_method %in% deprob_methods)

de_eval_measures_filtered <- de_eval_measures %>% dplyr::filter(is.na(mapping_type) | mapping_type == "de")
eval_measures_markers <- de_eval_measures_filtered %>% dplyr::filter(is.na(gene_set) | 
                                                                       gene_set == "markers")
eval_measures_full <- de_eval_measures_filtered %>% dplyr::filter(is.na(gene_set) | 
                                                                    gene_set == "full")
elist <- list(markers = eval_measures_markers, full = eval_measures_full)
de_plots <- lapply(elist, function(gs) {
  eval_measures_markers_de_melted <- gs %>% reshape2::melt(measure.vars = c("micro_f1",
                                                                            "accuracy"), 
                                                           variable.name = "measure", value.name = "value") %>% 
    dplyr::mutate_(.dots = setNames(list(lazyeval::interp(~factor(x), 
                                                          x = as.name("de_prob"))), "xval")) %>% dplyr::mutate(measure = plyr::mapvalues(measure, 
                                                                                                                                         c("micro_f1",
                                                                                                                                           "accuracy"),
                                                                                                                                         c("F1",
                                                                                                                                           "Accuracy")))
  method_ordering <- (eval_measures_markers_de_melted %>% 
                        dplyr::filter(de_prob == 0.35,
                                      measure == "F1") %>%
                        dplyr::group_by(clustering_method) %>%
                        dplyr::summarise(mean_value=mean(value, na.rm=TRUE)) %>%
                        dplyr::arrange(-mean_value))$clustering_method
  
  marker_de_plot <- ggplot(eval_measures_markers_de_melted %>%
                             dplyr::mutate(clustering_method=factor(clustering_method, levels = method_ordering)), 
                           aes(x = clustering_method, y = value, fill = clustering_method)) + 
    geom_boxplot(outlier.size = -1) + theme_bw() + theme_Publication() + 
    theme_nature() + stripped_theme() + facet_grid(measure~de_prob, 
                                                   scales = "free") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) + 
    xlab("Method") + ylab("Score") + 
    guides(fill = FALSE) 
  return(marker_de_plot)
})

de_plot_markers <- de_plots$markers + 
  guides(fill = FALSE) + 
  ggtitle("% of genes differentially expressed per cell type") + 
  scale_fill_manual(values = clust_methods_palette)

de_plot_full <- de_plots$full + 
  guides(fill = FALSE) + 
  ggtitle("% of genes differentially expressed per cell type") + 
  scale_fill_manual(values = clust_methods_palette)


de_plot_legend <- cellassign.utils::ggsimplelegend(names(clust_methods_palette),
                                                   colour_mapping = unname(clust_methods_palette),
                                                   legend_title = "Method", legend_rows = 2, fontsize = 7)
de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)

# Delta plots

delta_table <- de_deltas %>% 
  dplyr::filter(de_prob %in% delta_deprobs) %>%
  dplyr::mutate(de_prob = paste0("DE prob = ", de_prob))

rvals <- compute_pvals_subsets(delta_table,
                               facet_vars = c("de_prob", "clustering_method"),
                               formula = ~ true_delta + inferred_delta,
                               corfun = cor.test,
                               output = "estimate")

rval_labels <- rvals %>%
  dplyr::group_by(de_prob) %>%
  dplyr::summarise(r_label=as.character(as.expression(substitute(list(italic(R) == est1, italic(R[s]) == est2), list(est1 = format(estimate[clustering_method == "cellassign"], digits = 3),
                                                                                                                     est2 = format(estimate[clustering_method == "cellassign_shrinkage"], digits = 3))))))

delta_plots <- ggplot(delta_table, aes(x=true_delta, y=inferred_delta)) + 
  geom_point(aes(colour=clustering_method), alpha = 0.5) + 
  theme_Publication() + 
  theme_nature() + 
  stripped_theme() +
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(expand = c(0,0.05)) + 
  scale_y_continuous(expand = c(0,0.05)) + 
  xlab("True logFC") + 
  ylab("Inferred logFC") + 
  guides(colour = FALSE) + 
  facet_wrap(~ de_prob, ncol = length(delta_deprobs)) + 
  scale_colour_manual(values = clust_methods_palette) + 
  geom_text(data = rval_labels, aes(x=Inf, y=Inf, label=r_label), hjust = 1, vjust = 1, parse = TRUE,
            size = 0.35*8)

delta_plot_legend <- cellassign.utils::ggsimplelegend(unique(delta_table$clustering_method),
                                                      colour_mapping = unname(clust_methods_palette[unique(delta_table$clustering_method)]),
                                                      legend_title = "Method", legend_rows = 1, fontsize = 7)
delta_plot_legend <- cellassign.utils::extract_legend(delta_plot_legend)

# Final plot
de_plots_labeled <- cowplot::plot_grid(de_plot_full, de_plot_markers, de_plot_legend,
                                       labels = c('a', 'b', ''),
                                       ncol = 1,
                                       nrow = 3,
                                       rel_heights = c(1, 1, 0.2))

bottom_row <- cowplot::plot_grid(delta_plots,
                                 labels = c('c'),
                                 ncol = 1,
                                 nrow = 1)

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


