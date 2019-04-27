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
library(lmerTest)

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
parser$add_argument('--sce_fl_hgsc_merged', metavar='FILE', type = 'character',
                    help="SCE of merged FL and HGSC")
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
parser$add_argument('--fit_koh', metavar='FILE', type='character',
                    help="CellAssign fit to Koh data (full, no other)")
parser$add_argument('--fit_koh_hierarchical_top', metavar='FILE', type='character',
                    help="CellAssign fit to Koh data (hierarchical, top, no other)")
parser$add_argument('--fit_koh_hierarchical_bottom', metavar='FILE', type='character',
                    help="CellAssign fit to FL data (hierarchical, bottom, no other)")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_hgsc <- readRDS(args$sce_hgsc)
sce_fl <- readRDS(args$sce_follicular)
sce_follicular_hgsc_merged <- readRDS(args$sce_fl_hgsc_merged)
sce_follicular_hgsc_merged <- sce_follicular_hgsc_merged %>%
  scater::mutate(dataset = ifelse(str_detect(patient, "^VOA"), "HGSC", "FL"))

categorical_palettes <- cat_palettes()
combined_celltype_colours <- categorical_palettes$hgsc_fl_celltypes
cohort_colours <- categorical_palettes$cohort

koh_celltype_colours <- categorical_palettes$koh_celltype

fit_combined <- readRDS(args$fit_combined)
fit_hgsc <- readRDS(args$fit_hgsc)
fit_fl <- readRDS(args$fit_fl)

sce_follicular_lk <- readRDS(args$sce_follicular_lk)

fit_fl_noother <- readRDS(args$fit_fl_noother)
fit_fl_hierarchical_top <- readRDS(args$fit_fl_hierarchical_top)
fit_fl_hierarchical_bottom <- readRDS(args$fit_fl_hierarchical_bottom)

fit_koh <- readRDS(args$fit_koh)
fit_koh_hierarchical_top <- readRDS(args$fit_koh_hierarchical_top)
fit_koh_hierarchical_bottom <- readRDS(args$fit_koh_hierarchical_bottom)

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

create_hierarchical_probabilities <- function(fit_full, fit_top, fit_bottom,
                                              sce_full, sce_bottom,
                                              parent_type = "B cells",
                                              child_types = c("B cells (kappa)", "B cells (lambda)")) {
  gamma <- fit_bottom$mle_params$gamma
  
  full_info <- data.frame(Sample=sce_full$Sample,
                          Barcode=sce_full$Barcode,
                          celltype=fit_full$cell_type,
                          fit_full$mle_params$gamma,
                          check.names = FALSE)
  
  full_info_merged <- full_info
  full_info_merged[,parent_type] <- rowSums(full_info_merged[,child_types,drop=FALSE])
  full_info_merged <- full_info_merged %>%
    dplyr::mutate(celltype = plyr::mapvalues(celltype,
                                             from = child_types,
                                             to = rep(parent_type, length(child_types))))
  
  broad_info <- data.frame(Sample=sce_full$Sample,
                           Barcode=sce_full$Barcode,
                           celltype=fit_top$cell_type,
                           fit_top$mle_params$gamma,
                           check.names = FALSE)
  
  lk_info <- data.frame(Sample=sce_bottom$Sample,
                        Barcode=sce_bottom$Barcode,
                        celltype=fit_bottom$cell_type,
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
                       dplyr::filter(type == parent_type) %>%
                       dplyr::mutate(parent_prob=prob_hierarchical) %>%
                       dplyr::select(Sample, Barcode, parent_prob),
                     by = c("Sample", "Barcode")) %>%
    dplyr::mutate(prob_hierarchical=prob_hierarchical * parent_prob)
  
  
  broad_merged <- full_info_merged_melted %>%
    dplyr::inner_join(broad_info_melted)
  
  deep_merged <- full_info_melted %>%
    dplyr::inner_join(lk_info_melted)
  
  prob_df <- broad_merged %>%
    plyr::rbind.fill(deep_merged)
  
  return(prob_df)
}


fit_fl_hierarchical_bottom_renamed <- fit_fl_hierarchical_bottom
fit_fl_hierarchical_bottom_renamed$cell_type <- fit_fl_hierarchical_bottom_renamed$cell_type %>%
  plyr::mapvalues(from = c("IGKC", "IGLC"),
                  to = c("B cells (kappa)",
                         "B cells (lambda)"))
colnames(fit_fl_hierarchical_bottom_renamed$mle_params$gamma) <- colnames(fit_fl_hierarchical_bottom_renamed$mle_params$gamma) %>%
  plyr::mapvalues(from = c("IGKC", "IGLC"),
                  to = c("B cells (kappa)",
                         "B cells (lambda)"))

fl_hierarchical_probs <- create_hierarchical_probabilities(fit_fl_noother, 
                                                           fit_fl_hierarchical_top, 
                                                           fit_fl_hierarchical_bottom_renamed,
                                                           sce_fl, 
                                                           sce_follicular_lk,
                                                           parent_type = "B cells",
                                                           child_types = c("B cells (kappa)", "B cells (lambda)")) 

repeat_barcodes <- colData(sce_follicular_lk) %>%
  as.data.frame %>%
  dplyr::select(Sample, Barcode) %>%
  dplyr::mutate(repeat_barcode=1)

fl_hierarchical_probs <- fl_hierarchical_probs %>%
  dplyr::left_join(repeat_barcodes) %>%
  dplyr::filter(is.na(repeat_barcode) | (type != "B cells"))

mod <- lmer(prob_hierarchical ~ prob_full + (1|type), data = fl_hierarchical_probs)
pval <- unname(summary(mod)$coefficients[,5]['prob_full'])

eq <- substitute(italic(R)^2~"="~r2~", "~italic(P)~"="~p, 
                 list(r2 = format(with(fl_hierarchical_probs, cor(prob_hierarchical, prob_full)^2), digits = 3),
                      p = format(pval, digits = 3)))
eq <- as.character(as.expression(eq))

prob_plot_hierarchical_fl <- ggplot(fl_hierarchical_probs, aes(x=prob_hierarchical, y=prob_full)) + 
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


## Koh hierarchical plot

sce_koh$Sample <- "Temp"
sce_koh$Barcode <- paste0("Barcode", 1:ncol(sce_koh))

koh_hierarchical_probs <- create_hierarchical_probabilities(fit_koh, 
                                                            fit_koh_hierarchical_top, 
                                                            fit_koh_hierarchical_bottom,
                                                            sce_koh, 
                                                            sce_koh[,fit_koh_hierarchical_top$cell_type == "APS_MPS"],
                                                            parent_type = "APS_MPS",
                                                            child_types = c("APS", "MPS")) 

koh_hierarchical_probs <- koh_hierarchical_probs %>%
  dplyr::filter(type != "APS_MPS")

mod <- lmer(prob_hierarchical ~ prob_full + (1|type), data = koh_hierarchical_probs)
pval <- unname(summary(mod)$coefficients[,5]['prob_full'])

eq <- substitute(italic(R)^2~"="~r2~", "~italic(P)~"="~p, 
                 list(r2 = format(with(koh_hierarchical_probs, cor(prob_hierarchical, prob_full)^2), digits = 3),
                      p = format(pval, digits = 3)))
eq <- as.character(as.expression(eq))

prob_plot_hierarchical_koh <- ggplot(koh_hierarchical_probs, aes(x=prob_hierarchical, y=prob_full)) + 
  geom_point(aes(colour=type), alpha = 0.3) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_color_manual(values = koh_celltype_colours) + 
  xlab("Probability (hierarchical)") + 
  ylab("Probability (combined)") + 
  guides(colour = FALSE) + 
  annotate(geom = 'text', x = 0, y = Inf, hjust = -0.5, vjust = 1.5, label = eq, parse = TRUE,
           size = 3)

koh_celltype_legend <- cellassign.utils::ggsimplelegend(names(koh_celltype_colours),
                                                    colour_mapping = unname(koh_celltype_colours),
                                                    legend_title = "Celltype", legend_rows = 2, fontsize = 7)
koh_celltype_legend <- cellassign.utils::extract_legend(koh_celltype_legend)


## Plotting

hgsc_fl_legend_row <- cowplot::plot_grid(cohort_legend,
                                         celltype_legend,
                                         ncol = 2, 
                                         rel_widths = c(0.5, 0.5))

combined_hgsc_fl_plots <- cowplot::plot_grid(hgsc_fl_dimred_plots[[1]],
                                             hgsc_fl_dimred_plots[[2]],
                                             cohort_plot,
                                             prob_plot_hierarchical,
                                             labels = c('a', 'b', 'c', 'd'),
                                             nrow = 2,
                                             ncol = 2)

koh_combined_plot <- cowplot::plot_grid(prob_plot_hierarchical_koh,
                                        koh_celltype_legend,
                                        nrow = 2,
                                        rel_heights = c(0.85, 0.15))

bottom_row <- cowplot::plot_grid(koh_combined_plot,
                                 ncol = 2,
                                 rel_widths = c(0.5, 0.5),
                                 labels = c('e', ''))

final_plot <- cowplot::plot_grid(combined_hgsc_fl_plots,
                                 hgsc_fl_legend_row,
                                 bottom_row,
                                 nrow = 3,
                                 ncol = 1,
                                 rel_heights = c(0.9, 0.1, 0.5))


# Plot final plot
pdf(args$outfname, width = 10, height = 11.2, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


