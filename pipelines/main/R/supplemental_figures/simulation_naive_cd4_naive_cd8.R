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
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create supplemental wrong marker figure.")
parser$add_argument('--deprob_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for DE prob")
parser$add_argument('--delta_deprobs', type='double', nargs ='+',
                    help="DE probs to show delta plots for.")
parser$add_argument('--wrongmarker_result_dir', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker")
parser$add_argument('--wrongmarker_result_dir2', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker (with 20 genes per type)")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--method_description', type='character', metavar='FILE',
                    help="Method description")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir
delta_deprobs <- unlist(args$delta_deprobs)
deprob_methods <- unlist(args$deprob_methods)

wrongmarker_result_dir <- args$wrongmarker_result_dir
wrongmarker_result_dir2 <- args$wrongmarker_result_dir2
method_description <- read_yaml(args$method_description)
method_metadata <- plyr::rbind.fill(lapply(names(method_description), function(i) {
  data.frame(method_type=i, clustering_method=method_description[[i]])
}))

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$method_types[df_as_map(method_metadata, deprob_methods, from = "clustering_method", to = "method_type")]
names(clust_methods_palette) <- deprob_methods

# DE prob figure

deprob_result_files <- Sys.glob(file.path(deprob_result_dir, "evaluate", "*", "*", "*.tsv"))
de_eval_measures <- plyr::rbind.fill(lapply(deprob_result_files, function(f) {
  df <- fread(f)
  if ("gene_set" %in% colnames(df)) {
    return(df)
  } else {
    feature_type <- str_extract(f, "(markers|full)")
    return(data.frame(fread(f), gene_set=feature_type))
  }
}))
de_eval_measures <- de_eval_measures %>%
  dplyr::left_join(method_metadata)

delta_files <- Sys.glob(file.path(deprob_result_dir, "assign_celltypes_sce", "deltas", "*", "cellassign*.tsv"))
de_deltas <- plyr::rbind.fill(lapply(delta_files, function(f) {
  fread(f)
}))

## Only use methods that were selected
de_eval_measures <- de_eval_measures %>%
  dplyr::filter(clustering_method %in% deprob_methods)

## TODO: When reruns have been done with ARI and NMI, add those too
de_eval_measures_filtered <- de_eval_measures %>% 
  dplyr::filter(is.na(mapping_type) | mapping_type == "de",
                de_prob <= 0.55) 
eval_measures_markers <- de_eval_measures_filtered %>% dplyr::filter(is.na(gene_set) | 
                                                                       gene_set == "markers")
eval_measures_full <- de_eval_measures_filtered %>% dplyr::filter(is.na(gene_set) | 
                                                                    gene_set == "full")
elist <- list(markers = eval_measures_markers, full = eval_measures_full)

get_change_indexes <- function(x) {
  which(x[-1] != x[-length(x)])
}

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
  method_ordering_df <- eval_measures_markers_de_melted %>% 
    dplyr::filter(de_prob == 0.35,
                  measure == "F1") %>%
    dplyr::group_by(method_type, clustering_method) %>%
    dplyr::summarise(mean_value=mean(value, na.rm=TRUE)) %>%
    dplyr::arrange(method_type, -mean_value)
  
  method_ordering <- method_ordering_df$clustering_method
  
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
          wilcox_res <- wilcox.test(x_cast[,as.character(m1)], x_cast[,as.character(m2)], paired = TRUE)
          
          # Get direction of significance
          wilcox_res_greater <- wilcox.test(x_cast[,as.character(m1)], x_cast[,as.character(m2)], paired = TRUE, alternative = "greater")
          
          pval_max <- max(wilcox_res$p.value, wilcox_res_greater$p.value)
          return(data.frame(method1=m1, method2=m2, p.value=pval_max))
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
    geom_boxplot(outlier.size = -1, width = 0.6) + 
    geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() + 
    stripped_theme() + 
    facet_grid(measure~de_prob, 
               scales = "free") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.8),
          axis.title = element_text(size = 9.5, face = "bold"),
          plot.title = element_text(size = 11, face = "plain"),
          strip.text = element_text(size = 9.5, face = "plain")) + 
    xlab("Method") + ylab("Score") + 
    guides(fill = FALSE) + 
    geom_text(data = paired_pvals_cast,
              aes(label=cellassign, x=clustering_method, y=-0.22),
              colour = clust_methods_palette["cellassign"]) + 
    coord_cartesian(ylim = c(0, 1), clip = 'off') + 
    theme(panel.spacing.y = unit(2, "lines"))
  
  change_indexes <- get_change_indexes(method_ordering_df$method_type)
  for (x in change_indexes) {
    marker_de_plot <- marker_de_plot + 
      geom_vline(xintercept = x + 0.5, linetype = "dashed", alpha = 0.8)
  }
  
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


# Wrong marker figure (1)
wrongmarker_de_files <- Sys.glob(file.path(wrongmarker_result_dir, "evaluate", "celltypes", "*", "cellassign*.tsv"))
wm_eval_measures <- plyr::rbind.fill(lapply(wrongmarker_de_files, function(f) {
  fread(f)
}))

wm_plots <- plot_simulation_performance(wm_eval_measures,
                                        measures = c("micro_f1",
                                                     "accuracy"),
                                        display_measure_names = c("F1",
                                                                  "Accuracy"),
                                        x_var = "wrong_marker_proportion")

wm_plot_cellassign <- wm_plots$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = FALSE) + 
  scale_fill_manual(values = clust_methods_palette) + 
  theme(axis.title = element_text(size = 9.5, face = "bold"),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9.5, face = "plain"))

# Wrong marker figure (2)
# wm_eval_measures2 <- load_annotation_files(wrongmarker_result_dir2, pattern = "*_eval_measures.tsv")
# wm_delta_vals2 <- load_annotation_files(wrongmarker_result_dir2, pattern = "*_delta_compare.tsv")

wrongmarker_de_files2 <- Sys.glob(file.path(wrongmarker_result_dir2, "evaluate", "celltypes", "*", "cellassign*.tsv"))
wm_eval_measures2 <- plyr::rbind.fill(lapply(wrongmarker_de_files2, function(f) {
  fread(f)
}))

wm_plots2 <- plot_simulation_performance(wm_eval_measures2,
                                         measures = c("micro_f1",
                                                      "accuracy"),
                                         display_measure_names = c("F1",
                                                                   "Accuracy"),
                                         x_var = "wrong_marker_proportion")

wm_plot_cellassign2 <- wm_plots2$cellassign + 
  xlab("Proportion of incorrect entries in rho") + 
  guides(fill = FALSE) + 
  scale_fill_manual(values = clust_methods_palette) + 
  theme(axis.title = element_text(size = 9.5, face = "bold"),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9.5, face = "plain"))

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
  dplyr::summarise(r_label=as.character(as.expression(substitute(list(italic(R) == est1), list(est1 = format(estimate[clustering_method == "cellassign"], digits = 3))))))

delta_plots <- ggplot(delta_table, aes(x=true_delta, y=inferred_delta)) + 
  geom_point(alpha = 0.5) + 
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
  geom_text(data = rval_labels, aes(x=Inf, y=Inf, label=r_label), hjust = 1.05, vjust = 1.2, parse = TRUE,
            size = 0.35*9.5) + 
  theme(axis.title = element_text(size = 9.5, face = "bold"),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9.5, face = "plain"))


whitespace_width <- 0.4
rel_width <- 1 - whitespace_width
texts <- paste(c("p<0.001", "p<0.01", "p<0.05"), c("***", "**", "*"), sep = ":")
significance_legend <- gridExtra::arrangeGrob(grobs = lapply(texts, function(x) grid::textGrob(x, gp = grid::gpar(fontsize = 8))), layout_matrix = rbind(c(NA, 1:3, NA)),
                                              widths = c(whitespace_width, rep(rel_width/3, 3), whitespace_width))



de_plot_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$method_types),
                                                   colour_mapping = unname(categorical_palettes$method_types),
                                                   legend_title = "Method type", legend_rows = 1, fontsize = 7)
de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)

# Plot final plot

legend_row <- cowplot::plot_grid(de_plot_legend,
                                 significance_legend,
                                 ncol = 2, 
                                 rel_widths = c(0.5, 0.5))

wm_row <- cowplot::plot_grid(wm_plot_cellassign,
                             wm_plot_cellassign2,
                             ncol = 2,
                             rel_widths = c(0.5, 0.5), 
                             labels = c('c', 'd'))

final_plot <- cowplot::plot_grid(de_plot_markers,
                                 legend_row,
                                 delta_plots,
                                 wm_row,
                                 nrow = 4,
                                 labels = c('a', '', 'b', ''),
                                 rel_heights= c(1.25, 0.1, 1, 1))

pdf(args$outfname, width = 10, height = 11, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


