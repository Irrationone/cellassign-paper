# Simulation results for all, correlation-based mapping

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
parser$add_argument('--deprob_result_dir1', metavar='DIR', type='character',
                    help="Path to first simulation result directory for DE prob")
parser$add_argument('--deprob_result_dir2', metavar='DIR', type='character',
                    help="Path to second simulation result directory for DE prob")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--method_description', type='character', metavar='FILE',
                    help="Method description")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir1 <- args$deprob_result_dir1
deprob_result_dir2 <- args$deprob_result_dir2
deprob_methods <- unlist(args$deprob_methods)
method_description <- read_yaml(args$method_description)
method_metadata <- plyr::rbind.fill(lapply(names(method_description), function(i) {
  data.frame(method_type=i, clustering_method=method_description[[i]])
}))

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

clust_methods_palette <- categorical_palettes$method_types[df_as_map(method_metadata, deprob_methods, from = "clustering_method", to = "method_type")]
names(clust_methods_palette) <- deprob_methods

# DE prob figures
de_plots_combined <- lapply(c(deprob_result_dir1, deprob_result_dir2), function(deprob_result_dir) {
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
  
  ## Only use methods that were selected
  de_eval_measures <- de_eval_measures %>%
    dplyr::filter(clustering_method %in% deprob_methods)
  
  de_eval_measures_filtered <- de_eval_measures %>% dplyr::filter(is.na(mapping_type) | mapping_type == "correlation")
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
      geom_boxplot(outlier.size = -1) + 
      geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) + 
      theme_bw() + theme_Publication() + 
      theme_nature() + stripped_theme() + facet_grid(measure~de_prob, 
                                                     scales = "free") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.88)) + 
      xlab("Method") + ylab("Score") + 
      guides(fill = FALSE) + 
      geom_text(data = paired_pvals_cast,
                aes(label=cellassign, x=clustering_method, y=-0.15),
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
  
  de_plots$markers <- de_plots$markers + 
    guides(fill = FALSE) + 
    ggtitle("% of genes differentially expressed per cell type") + 
    scale_fill_manual(values = clust_methods_palette)
  
  de_plots$full <- de_plots$full + 
    guides(fill = FALSE) + 
    ggtitle("% of genes differentially expressed per cell type") + 
    scale_fill_manual(values = clust_methods_palette)
  
  return(de_plots)
})


de_plot_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$method_types),
                                                   colour_mapping = unname(categorical_palettes$method_types),
                                                   legend_title = "Method type", legend_rows = 1, fontsize = 7)
de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)

whitespace_width <- 0.3
rel_width <- 1 - whitespace_width
texts <- paste(c("p<0.001", "p<0.01", "p<0.05"), c("***", "**", "*"), sep = ":")
significance_legend <- gridExtra::arrangeGrob(grobs = lapply(texts, function(x) grid::textGrob(x, gp = grid::gpar(fontsize = 8))), layout_matrix = rbind(c(NA, 1:3, NA)),
                                              widths = c(whitespace_width, rep(rel_width/3, 3), whitespace_width))



# Final plot
legend_row <- cowplot::plot_grid(de_plot_legend,
                                 significance_legend,
                                 ncol = 2, 
                                 rel_widths = c(0.67, 0.33))

final_plot <- cowplot::plot_grid(de_plots_combined[[1]]$full, 
                                 de_plots_combined[[1]]$markers,
                                 de_plots_combined[[2]]$full,
                                 de_plots_combined[[2]]$markers, 
                                 legend_row,
                                 labels = c('a', 'b', 'c', 'd', ''),
                                 ncol = 1,
                                 nrow = 5,
                                 rel_heights = c(1, 1, 1, 1, 0.2))

# Plot final plot
pdf(args$outfname, width = 10, height = 15, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


