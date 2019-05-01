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
library(yaml)
library(ggbeeswarm)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create simulation figure.")
parser$add_argument('--deprob_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for DE prob")
parser$add_argument('--novel_extra_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for novel/extra celltype analysis")
parser$add_argument('--sce_liver', metavar='FILE', type='character',
                    help="SCE for liver")
parser$add_argument('--liver_marker_types', type='character', nargs = '+',
                    help="Celltypes with marker genes")
parser$add_argument('--cellassign_fit_novel8', metavar='FILE', type='character',
                    help="CellAssign novel liver fit")
parser$add_argument('--sce_mix_merged', metavar='FILE', type = 'character',
                    help="SCE of merged mixture data")
parser$add_argument('--fit_tian', metavar='FILE', type='character',
                    help="CellAssign fit to CellBench mixture")
parser$add_argument('--cell_lines', type='character', nargs = '+',
                    help="Cell lines to use")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--method_description', type='character', metavar='FILE',
                    help="Method description")
parser$add_argument('--dimreduce_type', type='character',
                    help="Dimreduce type", choices = c("PCA", "TSNE", "UMAP"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

deprob_result_dir <- args$deprob_result_dir
novel_extra_result_dir <- args$novel_extra_result_dir
cellassign_fit_novel8 <- readRDS(args$cellassign_fit_novel8)
liver_marker_types <- unlist(args$liver_marker_types)
sce_liver <- readRDS(args$sce_liver)

fit_tian <- readRDS(args$fit_tian)
cell_lines <- unlist(args$cell_lines)

sce_tian_mix <- readRDS(args$sce_mix_merged)

deprob_methods <- unlist(args$deprob_methods)
method_description <- read_yaml(args$method_description)
method_metadata <- plyr::rbind.fill(lapply(names(method_description), function(i) {
  data.frame(method_type=i, clustering_method=method_description[[i]])
}))

categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

liver_celltype_palette <- categorical_palettes$liver_celltypes
liver_celltype_palette <- c(liver_celltype_palette, 'Unassigned'='gray60')

clust_methods_palette <- categorical_palettes$method_types[df_as_map(method_metadata, deprob_methods, from = "clustering_method", to = "method_type")]
names(clust_methods_palette) <- deprob_methods

# DE prob figure
# de_eval_measures <- load_annotation_files(deprob_result_dir, pattern = "*_eval_measures.tsv")
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

#de_deltas <- load_annotation_files(deprob_result_dir, pattern = "*_delta_compare.tsv")
# delta_files <- Sys.glob(file.path(deprob_result_dir, "assign_celltypes_sce", "deltas", "*", "cellassign*.tsv"))
# de_deltas <- plyr::rbind.fill(lapply(delta_files, function(f) {
#   fread(f)
# }))

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

de_plot_markers <- de_plots$markers + 
  guides(fill = FALSE) + 
  ggtitle("% of genes differentially expressed per cell type") + 
  scale_fill_manual(values = clust_methods_palette)
  
de_plot_full <- de_plots$full + 
  guides(fill = FALSE) + 
  ggtitle("% of genes differentially expressed per cell type") + 
  scale_fill_manual(values = clust_methods_palette)


whitespace_width <- 0.4
rel_width <- 1 - whitespace_width
texts <- paste(c("p<0.001", "p<0.01", "p<0.05"), c("***", "**", "*"), sep = ":")
significance_legend <- gridExtra::arrangeGrob(grobs = lapply(texts, function(x) grid::textGrob(x, gp = grid::gpar(fontsize = 8))), layout_matrix = rbind(c(NA, 1:3, NA)),
                                              widths = c(whitespace_width, rep(rel_width/3, 3), whitespace_width))



de_plot_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$method_types),
                                                   colour_mapping = unname(categorical_palettes$method_types),
                                                   legend_title = "Method type", legend_rows = 1, fontsize = 7)
de_plot_legend <- cellassign.utils::extract_legend(de_plot_legend)


# Novel extra celltypes

novel_extra_result_files <- Sys.glob(file.path(novel_extra_result_dir, "evaluate", "celltypes", "*", "*.tsv"))
novel_extra_eval_measures <- plyr::rbind.fill(lapply(novel_extra_result_files, function(f) {
  df <- fread(f)
  if ("gene_set" %in% colnames(df)) {
    return(df)
  } else {
    feature_type <- str_extract(f, "(markers|full)")
    return(data.frame(fread(f), gene_set=feature_type))
  }
}))

novel_extra_eval_measures_melted <- novel_extra_eval_measures %>% 
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

novel_df <- novel_extra_eval_measures_melted %>%
  dplyr::filter(n_data_types >= n_marker_types) %>%
  dplyr::mutate(n_unknown_types=n_data_types - n_marker_types)

superset_df <- novel_extra_eval_measures_melted %>%
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


novel_min_val <- max(0, min(novel_df$value) - 0.02)
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
            aes(label=symbol, x=clustering_method, y=novel_min_val-0.15*(1-novel_min_val)),
            colour = clust_methods_palette["cellassign"]) + 
  coord_cartesian(ylim = c(novel_min_val, 1), clip = 'off') + 
  theme(panel.spacing.y = unit(2, "lines"),
        axis.title = element_text(size = 9.5, face = "bold"),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9.5, face = "plain")) + 
  scale_fill_manual(values = clust_methods_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Number of unknown cell types (out of 6)")

superset_min_val <- max(0, min(superset_df$value) - 0.02)
superset_plot <- ggplot(superset_df, 
                        aes(x = clustering_method, y = value, fill = clustering_method)) + 
  geom_boxplot(outlier.size = -1) + 
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) + 
  theme_bw() + theme_Publication() + 
  theme_nature() + stripped_theme() + facet_grid(measure~n_extra_types, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, vjust = 0.88)) + 
  xlab("Method") + ylab("Score") + 
  guides(fill = FALSE) + 
  coord_cartesian(ylim = c(superset_min_val, 1), clip = 'off') + 
  geom_text(data = superset_pvals %>% dplyr::rename(clustering_method=method2),
            aes(label=symbol, x=clustering_method, y=superset_min_val-0.15*(1-superset_min_val)),
            colour = clust_methods_palette["cellassign"]) +
  theme(panel.spacing.y = unit(2, "lines"),
        axis.title = element_text(size = 9.5, face = "bold"),
        plot.title = element_text(size = 11, face = "plain"),
        strip.text = element_text(size = 9.5, face = "plain")) + 
  scale_fill_manual(values = clust_methods_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Number of removed cell types in the data (out of 6)")


## Liver celltypes 

celltypes <- cellassign_fit_novel8$cell_type %>%
  plyr::mapvalues("other", "Unassigned")

ground_truth_celltypes <- sce_liver$celltype
ground_truth_celltypes[!ground_truth_celltypes %in% liver_marker_types] <- "Unassigned"

celltype_levels <- sort(unique(c(ground_truth_celltypes, celltypes)))
cont_table <- table(factor(ground_truth_celltypes, levels = celltype_levels), 
                    factor(celltypes, levels = celltype_levels))
acc <- sum(diag(cont_table))/sum(cont_table)
micro_f1 <- microF1(cont_table)
plot_label <- paste0("Acc = ",  format(acc, digits = 3), "\n", "F1 = ", format(micro_f1, digits = 3))
  
sce_liver$assigned_cluster <- celltypes

liver_novel_plot <- plotReducedDim(sce_liver,
                                   use_dimred = args$dimreduce_type,
                                   colour_by = "assigned_cluster",
                                   point_alpha = 0.4, 
                                   point_size = 1.5)
liver_novel_plot <- liver_novel_plot + 
  guides(colour = FALSE,
         shape = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = liver_celltype_palette) + 
  guides(fill = FALSE) + 
  ggtitle("CellAssign (reduced markers)") + 
  annotate(geom = 'text', x = Inf, y = Inf, hjust = 1, vjust = 1.5, label = plot_label, parse = FALSE,
           size = 2.5)

master_label_plot <- plotReducedDim(sce_liver,
                                    use_dimred = args$dimreduce_type,
                                    colour_by = "celltype",
                                    point_alpha = 0.4, 
                                    point_size = 1.5) +
  guides(colour = FALSE,
         shape = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = liver_celltype_palette) + 
  guides(fill = FALSE) + 
  ggtitle("Cell type")

celltype_legend <- cellassign.utils::ggsimplelegend(names(liver_celltype_palette),
                                                    colour_mapping = unname(liver_celltype_palette),
                                                    legend_title = "Celltype", legend_rows = 3, fontsize = 7)
celltype_legend <- cellassign.utils::extract_legend(celltype_legend)




## CellBench data

for (cl in cell_lines) {
  colData(sce_tian_mix)[,cl][sce_tian_mix$cell_line == cl] <- 9
  colData(sce_tian_mix)[,cl][sce_tian_mix$cell_line != cl] <- 0
}

cellassign_probs <- reshape2::melt(fit_tian$mle_params$gamma) %>%
  dplyr::rename(cell_id=Var1, cell_line=Var2, cellassign_prob=value)

true_probs <- colData(sce_tian_mix)[,c("sample_id", cell_lines)] %>%
  as.data.frame %>%
  dplyr::mutate(cell_id=1:n()) %>%
  reshape2::melt(id.vars = c("cell_id", "sample_id"), variable.name = "cell_line", value.name = "num_cells")

valid_cells <- true_probs %>%
  dplyr::group_by(sample_id, cell_id) %>%
  dplyr::summarise(total_cells=sum(num_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(total_cells == "9")

prob_df <- true_probs %>%
  dplyr::left_join(cellassign_probs) %>%
  dplyr::filter(cell_id %in% valid_cells$cell_id,
                sample_id == "mixture") %>%
  dplyr::mutate(clustering_method = "cellassign")

## CellBench probability plot
cellbench_prob_plot <- ggplot(prob_df, aes(x=factor(num_cells), y=cellassign_prob)) + 
  geom_boxplot(width = 0.5, outlier.size = -1, aes(fill = clustering_method)) + 
  geom_quasirandom(alpha = 0.2, width = 0.25) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  facet_wrap(~ cell_line, ncol = 1) + 
  stripped_theme(strip_face = "bold") + 
  guides(fill = FALSE) + 
  scale_fill_manual(values = clust_methods_palette) + 
  xlab("# cells in 9-cell barcode") + 
  ylab("CellAssign probability") + 
  theme(plot.margin = unit(c(5, 5, 5, 5), "mm"))



# Final plot
legend_row <- cowplot::plot_grid(de_plot_legend,
                                 significance_legend,
                                 ncol = 2, 
                                 rel_widths = c(0.5, 0.5))

novel_superset_row <- cowplot::plot_grid(superset_plot,
                                         novel_plot,
                                         ncol = 2,
                                         labels = c('b', 'c'))

de_plot_rows <- cowplot::plot_grid(de_plot_full, #de_plot_legend,
                                 novel_superset_row,
                                 legend_row,
                                 labels = c('a', '', ''),
                                 ncol = 1,
                                 nrow = 3,
                                 rel_heights = c(1, 1, 0.1))


liver_plots <- cowplot::plot_grid(master_label_plot,
                                  liver_novel_plot,
                                  ncol = 2, 
                                  labels = c('d', 'e'))

liver_labeled <- cowplot::plot_grid(liver_plots,
                                    celltype_legend,
                                    nrow = 2, 
                                    rel_heights = c(0.75, 0.25))

lower_row <- cowplot::plot_grid(liver_labeled,
                                cellbench_prob_plot,
                                labels = c('', 'f'),
                                ncol = 2, 
                                rel_widths = c(0.5, 0.5))

final_plot <- cowplot::plot_grid(de_plot_rows,
                                 lower_row,
                                 nrow = 2,
                                 rel_heights = c(0.7, 0.3))

# Plot final plot
pdf(args$outfname, width = 10, height = 13, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


