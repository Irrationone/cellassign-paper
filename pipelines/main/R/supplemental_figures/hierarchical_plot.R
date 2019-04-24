# Hierarchical plot

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
library(ggrastr)
library(mclust)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Hierarchical plot")
parser$add_argument('--sce_hgsc', metavar='FILE', type='character',
                    help="SCE of HGSC")
parser$add_argument('--sce_follicular', metavar='FILE', type = 'character',
                    help="SCE of FL")
parser$add_argument('--sce_follicular_lk', metavar='FILE', type = 'character',
                    help="SCE of lambda-kappa on nonmalignant")
parser$add_argument('--fit_combined', metavar='FILE', type='character',
                    help="CellAssign fit to combined data")
parser$add_argument('--fit_hgsc', metavar='FILE', type='character',
                    help="CellAssign fit to HGSC data")
parser$add_argument('--fit_fl', metavar='FILE', type='character',
                    help="CellAssign fit to FL data, full")
parser$add_argument('--fit_fl_noother', metavar='FILE', type='character',
                    help="CellAssign fit to FL data (full, no other)")
parser$add_argument('--fit_fl_hierarchical_top', metavar='FILE', type='character',
                    help="CellAssign fit to FL data (hierarchical, top, no other)")
parser$add_argument('--fit_fl_hierarchical_bottom', metavar='FILE', type='character',
                    help="CellAssign fit to FL data (hierarchical, bottom, no other)")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_hgsc <- readRDS(args$sce_hgsc)
sce_fl <- readRDS(args$sce_fl)

categorical_palettes <- cat_palettes()
combined_celltype_colours <- categorical_palettes$hgsc_fl_celltypes
cohort_colours <- categorical_palettes$cohort

fit_combined <- readRDS(args$fit_combined)
fit_hgsc <- readRDS(args$fit_hgsc)
fit_fl <- readRDS(args$fit_fl)

sce_follicular_lk <- readRDS(args$sce_follicular_lk)

fit_fl_noother <- readRDS(args$fit_fl_noother)
fit_fl_hierarchical_top <- readRDS(args$fit_fl_hierarchical_top)
fit_fl_hierarchical_bottom <- readRDS(args$fit_fl_hierarchical_bottom)

merge_sces <- function(sces) {
  coldata_union <- Reduce(union, lapply(sces, function(x) {
    colnames(colData(x))
  }))
  
  common_genes <- Reduce(intersect, lapply(sces, function(x) {
    rownames(x)
  }))
  
  sces <- lapply(sces, function(x) {
    missing_cols <- setdiff(coldata_union, colnames(colData(x)))
    rowdat_cols <- setdiff(colnames(rowData(x)), c("mean_counts", "log10_mean_counts",
                                                   "n_cells_by_counts", "pct_dropout_by_counts", "total_counts",
                                                   "log10_total_counts"))
    rowData(x) <- rowData(x)[,rowdat_cols]
    colData(x)[,missing_cols] <- NA
    colData(x) <- colData(x)[,coldata_union]
    x <- x[common_genes,]
    return(x)
  })
  
  sce <- Reduce(cbind, sces)
  
  # Recompute size factors
  qclust <- quickCluster(sce, min.size = 30)
  sce <- computeSumFactors(sce, clusters = qclust)
  
  sce$size_factor <- sizeFactors(sce)
  return(sce)
}

sce_follicular_hgsc_merged <- merge_sces(list(sce_fl, sce_hgsc))
sce_follicular_hgsc_merged <- normalize(sce_follicular_hgsc_merged)
sce_follicular_hgsc_merged <- runPCA(sce_follicular_hgsc_merged, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
sce_follicular_hgsc_merged <- runTSNE(sce_follicular_hgsc_merged, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
sce_follicular_hgsc_merged <- runUMAP(sce_follicular_hgsc_merged, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)

sce_follicular_hgsc_merged <- sce_follicular_hgsc_merged %>%
  scater::mutate(dataset = ifelse(str_detect(patient, "^VOA"), "HGSC", "FL"))

sce_follicular_hgsc_merged$celltype_combined <- fit_combined$cell_type %>%
  plyr::mapvalues(from = c("B cells (lambda)", "B cells (kappa)", "other"), to = c("B cells", "B cells", "Unassigned"))
sce_follicular_hgsc_merged$celltype_individual <- c(fit_fl$cell_type, fit_hgsc$cell_type) %>%
  plyr::mapvalues(from = c("B cells (lambda)", "B cells (kappa)", "other"), to = c("B cells", "B cells", "Unassigned"))

plot_titles <- c("CellAssign (separate)", "CellAssign (combined)")
categories <- c("celltype_individual", "celltype_combined")

hgsc_fl_dimred_plots <- lapply(seq_along(categories), function(i) {
  x <- categories[i]
  p <- plotReducedDim(sce_follicular_hgsc_merged,
                      use_dimred = args$dimreduce_type,
                      colour_by = x,
                      point_alpha = 0.4, 
                      point_size = 1.5) +
    guides(colour = FALSE,
           shape = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() + 
    scale_fill_manual(values = combined_celltype_colours) + 
    guides(fill = FALSE) + 
    ggtitle(plot_titles[i])
  return(p)
})

percent_equal <- sum(sce_follicular_hgsc_merged$celltype_combined == sce_follicular_hgsc_merged$celltype_individual)/ncol(sce_follicular_hgsc_merged)
ari <- adjustedRandIndex(sce_follicular_hgsc_merged$celltype_combined, sce_follicular_hgsc_merged$celltype_individual)

plot_label <- paste0("% equal = ",  format(percent_equal, digits = 3), "\n", "ARI = ", format(ari, digits = 3))
hgsc_fl_dimred_plots[[2]] <- hgsc_fl_dimred_plots[[2]] + annotate(geom = 'text', x = Inf, y = Inf, hjust = 1, vjust = 1.5, label = plot_label, parse = FALSE,
                                                                  size = 2.5)


cohort_plot <- plotReducedDim(sce_follicular_hgsc_merged,
                              use_dimred = args$dimreduce_type,
                              colour_by = "dataset",
                              point_alpha = 0.4, 
                              point_size = 1.5) +
  guides(colour = FALSE,
         shape = FALSE) + 
  scale_fill_manual(values = cohort_colours) +
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  guides(fill = FALSE) + 
  ggtitle("Cohort")

## FL hierarchical 

gamma <- fit_fl_hierarchical_bottom$mle_params$gamma
colnames(gamma) <- colnames(gamma) %>%
  plyr::mapvalues(from = c("IGKC", "IGLC"),
                  to = c("B cells (kappa)",
                         "B cells (lambda)"))

full_info <- data.frame(Sample=sce_fl$Sample,
                        Barcode=sce_fl$Barcode,
                        celltype=fit_fl_noother$cell_type,
                        fit_fl_noother$mle_params$gamma,
                        check.names = FALSE)

full_info_merged <- full_info %>%
  dplyr::mutate(`B cells` = `B cells (lambda)` + `B cells (kappa)`,
                celltype = plyr::mapvalues(celltype,
                                           from = c("B cells (kappa)", "B cells (lambda)"),
                                           to = c("B cells", "B cells")))

broad_info <- data.frame(Sample=sce_fl$Sample,
                         Barcode=sce_fl$Barcode,
                         celltype=fit_fl_hierarchical_top$cell_type,
                         fit_fl_hierarchical_top$mle_params$gamma,
                         check.names = FALSE)

lk_info <- data.frame(Sample=sce_follicular_lk$Sample,
                      Barcode=sce_follicular_lk$Barcode,
                      celltype=plyr::mapvalues(fit_fl_hierarchical_bottom$cell_type, 
                                               from = c("IGKC", "IGLC"),
                                               to = c("B cells (kappa)",
                                                      "B cells (lambda)")),
                      gamma,
                      check.names = FALSE)


full_info_melted <- reshape2::melt(full_info, 
                                   id.vars = c("Sample", "Barcode", "celltype"),
                                   variable.name = "type", 
                                   value.name = "prob_full")
full_info_merged_melted <- reshape2::melt(full_info_merged, 
                                          id.vars = c("Sample", "Barcode", "celltype"),
                                          variable.name = "type", 
                                          value.name = "prob_full")
broad_info_melted <- reshape2::melt(broad_info, 
                                    id.vars = c("Sample", "Barcode", "celltype"),
                                    variable.name = "type", 
                                    value.name = "prob_hierarchical")
lk_info_melted <- reshape2::melt(lk_info, 
                                 id.vars = c("Sample", "Barcode", "celltype"),
                                 variable.name = "type", 
                                 value.name = "prob_hierarchical") %>%
  dplyr::left_join(broad_info_melted %>% 
                     dplyr::filter(type == "B cells") %>%
                     dplyr::mutate(b_prob=prob_hierarchical) %>%
                     dplyr::select(Sample, Barcode, b_prob),
                   by = c("Sample", "Barcode")) %>%
  dplyr::mutate(prob_hierarchical=prob_hierarchical * b_prob)



broad_merged <- full_info_merged_melted %>%
  dplyr::inner_join(broad_info_melted)

deep_merged <- full_info_melted %>%
  dplyr::inner_join(lk_info_melted)

prob_df <- broad_merged %>%
  plyr::rbind.fill(deep_merged)


eq <- substitute(italic(R)^2~"="~r2, 
                 list(r2 = format(with(prob_df, cor(prob_hierarchical, prob_full)^2), digits = 3)))
eq <- as.character(as.expression(eq))

prob_plot_hierarchical <- ggplot(prob_df, aes(x=prob_hierarchical, y=prob_full)) + 
  geom_point(aes(colour=type), alpha = 0.3) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_color_manual(values = combined_celltype_colours) + 
  xlab("Probability (hierarchical)") + 
  ylab("Probability (combined)") + 
  guides(colour = FALSE) + 
  annotate(geom = 'text', x = 0, y = Inf, hjust = -0.5, vjust = 1.5, label = eq, parse = TRUE,
           size = 3)


celltype_legend <- cellassign.utils::ggsimplelegend(names(combined_celltype_colours),
                                                    colour_mapping = unname(combined_celltype_colours),
                                                    legend_title = "Celltype", legend_rows = 4, fontsize = 7)
celltype_legend <- cellassign.utils::extract_legend(celltype_legend)

cohort_legend <- cellassign.utils::ggsimplelegend(names(cohort_colours),
                                                  colour_mapping = unname(cohort_colours),
                                                  legend_title = "Cohort", legend_rows = 1, fontsize = 7)
cohort_legend <- cellassign.utils::extract_legend(cohort_legend)

legend_row <- cowplot::plot_grid(cohort_legend,
                                 celltype_legend,
                                 ncol = 2, 
                                 rel_widths = c(0.4, 0.6))

combined_plots <- cowplot::plot_grid(hgsc_fl_dimred_plots[[1]],
                                     hgsc_fl_dimred_plots[[2]],
                                     cohort_plot,
                                     prob_plot_hierarchical,
                                     nrow = 2,
                                     ncol = 2)

final_plot <- cowplot::plot_grid(combined_plots,
                                 legend_row,
                                 nrow = 2,
                                 rel_heights = c(0.9, 0.1))


# Plot final plot
pdf(args$outfname, width = 10, height = 11.2, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


