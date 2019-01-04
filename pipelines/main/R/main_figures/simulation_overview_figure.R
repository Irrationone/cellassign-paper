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
parser$add_argument('--delta_deprobs', type='double', nargs ='+',
                    help="DE probs to show delta plots for.")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir
wrongmarker_result_dir <- args$wrongmarker_result_dir
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

## TODO: When reruns have been done with ARI and NMI, add those too
de_eval_measures_filtered <- de_eval_measures %>% 
  dplyr::filter(is.na(mapping_type) | mapping_type == "de",
                de_prob <= 0.45) # Just do up to 0.45
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
  
  df <- eval_measures_markers_de_melted %>%
    dplyr::mutate(clustering_method=factor(clustering_method, levels = method_ordering))
  
  paired_pvals <- df %>%
    plyr::ddply(plyr:::.(measure, de_prob), function(x) {
      x_cast <- reshape2::dcast(x, 
                                formula = seed ~ clustering_method, 
                                value.var = "value")
      clustering_methods <- df$clustering_method %>%
        unique
      cellassign_methods <- clustering_methods[str_detect(clustering_methods, "cellassign")]
      other_methods <- setdiff(clustering_methods, cellassign_methods)
      
      pvals <- plyr::rbind.fill(lapply(cellassign_methods, function(m1) {
        res <- plyr::rbind.fill(lapply(other_methods, function(m2) {
          wilcox_res <- wilcox.test(x_cast[,m1], x_cast[,m2], paired = TRUE)
          return(data.frame(method1=m1, method2=m2, p.value=wilcox_res$p.value))
        }))
        return(res)
      }))
      return(pvals)
    })
  
  paired_pvals <- paired_pvals %>%
    dplyr::mutate(p.adjust = p.adjust(paired_pvals$p.value, method = 'fdr'),
                  symbol=c("***", "**", "*", "")[.bincode(p.adjust, c(0, 1e-3, 1e-2, 0.05, 1))])
  
  df_extreme_vals <- df %>%
    dplyr::group_by(de_prob, clustering_method, measure) %>% 
    dplyr::summarise(max_val=max(value, na.rm=TRUE))
  
  paired_pvals_cast <- paired_pvals %>%
    reshape2::dcast(formula = measure + de_prob + method2 ~ method1, 
                    value.var = "symbol", 
                    fun.aggregate = function(x) x[1]) %>%
    dplyr::rename(clustering_method=method2) %>%
    dplyr::left_join(df_extreme_vals)
  
  marker_de_plot <- ggplot(df, 
                           aes(x = clustering_method, y = value, fill = clustering_method)) + 
    geom_boxplot(outlier.size = 0.4) + theme_bw() + theme_Publication() + 
    theme_nature() + stripped_theme() + facet_grid(measure~de_prob, 
                                                   scales = "free") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.8)) + 
    xlab("Method") + ylab("Score") + 
    guides(fill = FALSE) + 
    geom_text(data = paired_pvals_cast,
              aes(label=cellassign, x=clustering_method, y=-0.15),
              colour = clust_methods_palette["cellassign"]) + 
    geom_text(data = paired_pvals_cast,
              aes(label=cellassign_shrinkage, x=clustering_method, y=-0.22),
              colour = clust_methods_palette["cellassign_shrinkage"]) + 
    coord_cartesian(ylim = c(0, 1), clip = 'off') + 
    theme(panel.spacing.y = unit(2, "lines"))
    
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


whitespace_width <- 0.7
rel_width <- 1 - whitespace_width
texts <- paste(c("p<0.001", "p<0.01", "p<0.05"), c("***", "**", "*"), sep = ":")
significance_legend <- gridExtra::arrangeGrob(grobs = lapply(texts, function(x) grid::textGrob(x, gp = grid::gpar(fontsize = 8))), layout_matrix = rbind(c(NA, 1:3, NA)),
                                              widths = c(whitespace_width, rep(rel_width/3, 3), whitespace_width))



# 
# 
# de_plot_legend <- cellassign.utils::ggsimplelegend(names(clust_methods_palette),
#                                                    colour_mapping = unname(clust_methods_palette),
#                                                    legend_title = "Method", legend_rows = 2, fontsize = 7)
# de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)

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
  theme_bw() +
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
  geom_text(data = rval_labels, aes(x=Inf, y=Inf, label=r_label), hjust = 1.05, vjust = 1.2, parse = TRUE,
            size = 0.35*8)

delta_plot_legend <- cellassign.utils::ggsimplelegend(unique(delta_table$clustering_method),
                                                      colour_mapping = unname(clust_methods_palette[unique(delta_table$clustering_method)]),
                                                      legend_title = "Method", legend_rows = 1, fontsize = 7)
delta_plot_legend <- cellassign.utils::extract_legend(delta_plot_legend)
  

# Wrong marker figure
wm_eval_measures <- load_annotation_files(wrongmarker_result_dir, pattern = "*_eval_measures.tsv")
wm_delta_vals <- load_annotation_files(wrongmarker_result_dir, pattern = "*_delta_compare.tsv")

wm_plots <- plot_simulation_performance(wm_eval_measures %>%
                                          dplyr::filter(max_genes == 5), 
                                        measures = c("micro_f1",
                                                     "accuracy"),
                                        display_measure_names = c("F1",
                                                                  "Accuracy"),
                                        x_var = "wrong_marker_proportion")

wm_plot_cellassign <- wm_plots$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = FALSE) + 
  scale_fill_manual(values = clust_methods_palette)


# Final plot
de_plots_labeled <- cowplot::plot_grid(de_plot_full, de_plot_markers, #de_plot_legend,
                                       significance_legend,
                                       labels = c('a', 'b', ''),
                                       ncol = 1,
                                       nrow = 3,
                                       rel_heights = c(1, 1, 0.05))

bottom_row <- cowplot::plot_grid(delta_plots, wm_plot_cellassign,
                                 labels = c('c', 'd'),
                                 ncol = 2,
                                 nrow = 1,
                                 rel_widths = c(0.5, 0.5))

final_plot <- cowplot::plot_grid(de_plots_labeled, 
                                 bottom_row,
                                 delta_plot_legend,
                                 labels = c('', '', ''), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(0.75, 0.25, 0.05))


# Plot final plot
pdf(args$outfname, width = 10, height = 13, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


