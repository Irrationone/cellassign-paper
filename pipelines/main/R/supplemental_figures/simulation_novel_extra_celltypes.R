# Novel cell type simulation plots

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
                    help="Path to first simulation result directory for DE prob")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$clustering_methods

# DE prob figures
deprob_result_files <- Sys.glob(file.path(deprob_result_dir, "evaluate", "celltypes", "*", "*.tsv"))
de_eval_measures <- plyr::rbind.fill(lapply(deprob_result_files, function(f) {
  df <- fread(f)
  if ("gene_set" %in% colnames(df)) {
    return(df)
  } else {
    feature_type <- str_extract(f, "(markers|full)")
    return(data.frame(fread(f), gene_set=feature_type))
  }
}))

de_eval_measures_melted <- de_eval_measures %>% 
  reshape2::melt(measure.vars = c("micro_f1",
                                  "accuracy"), 
                 variable.name = "measure", value.name = "value") %>% 
  dplyr::mutate_(.dots = setNames(list(lazyeval::interp(~factor(x), 
                                                        x = as.name("de_prob"))), "xval")) %>% 
  dplyr::mutate(measure = plyr::mapvalues(measure, 
                                          c("micro_f1",
                                            "accuracy"),
                                          c("F1",
                                            "Accuracy")))

novel_df <- de_eval_measures_melted %>%
  dplyr::filter(n_data_types >= n_marker_types) %>%
  dplyr::mutate(n_unknown_types=n_data_types - n_marker_types)

superset_df <- de_eval_measures_melted %>%
  dplyr::filter(n_data_types <= n_marker_types) %>%
  dplyr::mutate(n_extra_types=n_marker_types - n_data_types)

compute_pvals <- function(x) {
  x_cast <- reshape2::dcast(x, 
                            formula = seed ~ clustering_method, 
                            value.var = "value")
  clustering_methods <- setdiff(colnames(x_cast), "seed")
  cellassign_methods <- clustering_methods[str_detect(clustering_methods, "cellassign")]
  other_methods <- setdiff(clustering_methods, cellassign_methods)
  
  pvals <- plyr::rbind.fill(lapply(cellassign_methods, function(m1) {
    res <- plyr::rbind.fill(lapply(other_methods, function(m2) {
      wilcox_res <- wilcox.test(x_cast[,as.character(m1)], x_cast[,as.character(m2)], paired = TRUE)
      
      # Get direction of significance
      wilcox_res_greater <- wilcox.test(x_cast[,as.character(m1)], x_cast[,as.character(m2)], paired = TRUE, alternative = "greater")
      
      pval_max <- max(wilcox_res$p.value, wilcox_res_greater$p.value)
      return(data.frame(method1=m1, method2=m2, p.value=pval_max))
    }))
    return(res)
  }))
  return(pvals)
}

novel_pvals <- novel_df %>%
  plyr::ddply(plyr:::.(measure, n_unknown_types), function(x) compute_pvals(x)) %>%
  dplyr::mutate(p.adjust = p.adjust(p.value, method = 'fdr'),
                symbol=c("***", "**", "*", "")[.bincode(p.adjust, c(0, 1e-3, 1e-2, 0.05, 1))])

superset_pvals <- superset_df %>%
  plyr::ddply(plyr:::.(measure, n_extra_types), function(x) compute_pvals(x)) %>%
  dplyr::mutate(p.adjust = p.adjust(p.value, method = 'fdr'),
                symbol=c("***", "**", "*", "")[.bincode(p.adjust, c(0, 1e-3, 1e-2, 0.05, 1))])

superset_df_extreme <- superset_df %>%
  dplyr::group_by(measure) %>%
  dplyr::summarise(min_measure=min(value)) %>%
  dplyr::ungroup()

novel_plot <- ggplot(novel_df, 
                         aes(x = clustering_method, y = value, fill = clustering_method)) + 
  geom_boxplot(outlier.size = -1) + 
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) + 
  theme_bw() + theme_Publication() + 
  theme_nature() + stripped_theme() + facet_grid(measure~n_unknown_types) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.88)) + 
  xlab("Method") + ylab("Score") + 
  guides(fill = FALSE) + 
  geom_text(data = novel_pvals %>% dplyr::rename(clustering_method=method2),
            aes(label=symbol, x=clustering_method, y=-0.15),
            colour = clust_methods_palette["cellassign"]) + 
  coord_cartesian(ylim = c(0, 1), clip = 'off') + 
  theme(panel.spacing.y = unit(2, "lines")) + 
  scale_fill_manual(values = clust_methods_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Number of unknown cell types (out of 6)")


superset_plot <- ggplot(superset_df, 
                     aes(x = clustering_method, y = value, fill = clustering_method)) + 
  geom_boxplot(outlier.size = -1) + 
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) + 
  theme_bw() + theme_Publication() + 
  theme_nature() + stripped_theme() + facet_grid(measure~n_extra_types, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.88)) + 
  xlab("Method") + ylab("Score") + 
  guides(fill = FALSE) + 
  coord_cartesian(ylim = c(0, 1), clip = 'off') + 
  geom_text(data = superset_pvals %>% dplyr::rename(clustering_method=method2) %>% dplyr::left_join(superset_df_extreme),
            aes(label=symbol, x=clustering_method, y=-0.15),
            colour = clust_methods_palette["cellassign"]) +
  theme(panel.spacing.y = unit(2, "lines")) + 
  scale_fill_manual(values = clust_methods_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Number of removed cell types in the data (out of 6)")

whitespace_width <- 0.7
rel_width <- 1 - whitespace_width
texts <- paste(c("p<0.001", "p<0.01", "p<0.05"), c("***", "**", "*"), sep = ":")
significance_legend <- gridExtra::arrangeGrob(grobs = lapply(texts, function(x) grid::textGrob(x, gp = grid::gpar(fontsize = 8))), layout_matrix = rbind(c(NA, 1:3, NA)),
                                              widths = c(whitespace_width, rep(rel_width/3, 3), whitespace_width))


# Final plot
final_plot <- cowplot::plot_grid(novel_plot,
                                 superset_plot,
                                 significance_legend,
                                 labels = c('a', 'b', ''),
                                 ncol = 1,
                                 nrow = 3,
                                 rel_heights = c(1, 1, 0.2))

# Plot final plot
pdf(args$outfname, width = 10, height = 13, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


