# Figure looking in-depth at nonmalignant cells in FL

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

parser <- ArgumentParser(description = "Create nonmalignant figure for FL")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_tcell', metavar='FILE', type='character',
                    help="Path to T cell-filtered SCE")
parser$add_argument('--sce_bcell', metavar='FILE', type='character',
                    help="Path to B cell-filtered SCE")
parser$add_argument('--sce_scvis_merged', metavar='FILE', type='character',
                    help="Path to scvis-merged SCE")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--tcell_labels', type='character', nargs='+',
                    help="Cell type labels of T cells")
parser$add_argument('--bcell_labels', type='character', nargs='+',
                    help="Cell type labels of B cells")
parser$add_argument('--azizi_signatures', metavar = 'FILE', type='character',
                    help="XLS file of Azizi et al. Table S4")
parser$add_argument('--de_cytotoxic_dir', type='character',
                    help="DE results for nonmalignant B cells")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
sce_tcell <- readRDS(args$sce_tcell)
sce_bcell <- readRDS(args$sce_bcell)
sce_scvis_merged <- readRDS(args$sce_scvis_merged)
de_cytotoxic_dir <- args$de_cytotoxic_dir

tcell_labels <- unlist(args$tcell_labels)
bcell_labels <- unlist(args$bcell_labels)

broad_celltype_map <- c(rep("T cell", length(tcell_labels)), 
                        rep("B cell", length(bcell_labels)),
                        "Unassigned")
names(broad_celltype_map) <- c(tcell_labels, bcell_labels, "Unassigned")


categorical_palettes <- cat_palettes()
heatmap_heat_colours <- heat_colour_gradient()
gradient_colours <- scrna_expression_gradient()

# Process SCEs
sce_tcell <- sce_tcell %>%
  scater::filter(celltype_full %in% c(tcell_labels, "other")) %>%
  scater::mutate(celltype_full=factor(plyr::mapvalues(celltype_full,
                                                      "other",
                                                      "Unassigned")))


sce_bcell <- sce_bcell %>%
  scater::filter(celltype_full %in% c(bcell_labels, "other")) %>%
  scater::mutate(celltype_full=factor(plyr::mapvalues(celltype_full,
                                                      "other",
                                                      "Unassigned")))

sce_scvis_merged <- sce_scvis_merged %>%
  scater::mutate(celltype_full=factor(plyr::mapvalues(celltype_full,
                                                      "other",
                                                      "Unassigned"))) %>%
  scater::mutate(category=ifelse(set_name == "follicular",
                                broad_celltype_map[as.character(celltype_full)],
                                 "RLN"),
                 category_full=ifelse(category == "B cell" & malignant_status_manual == "malignant",
                                      "B cell (malignant)",
                                      category)) 


# B cell dr figure
dr_bcell <- plotReducedDim(sce_bcell, 
                           use_dimred = "UMAP", 
                           colour_by = "celltype_full",
                           point_alpha = 0.2, 
                           add_ticks = FALSE,
                           point_size = 0.75)
dr_bcell$layers[[1]]$aes_params$colour <- NULL
dr_bcell$layers[[1]]$aes_params$shape <- 16
dr_bcell$layers[[1]]$mapping$colour <- dr_bcell$layers[[1]]$mapping$fill

dr_bcell <- dr_bcell + 
  guides(colour = guide_legend(title = "Celltype", override.aes = list(alpha = 1,
                                                                       size = 2), ncol = 2),
         fill = FALSE,
         shape = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$celltype)


dr_bcell_sample <- plotReducedDim(sce_bcell, 
                                  use_dimred = "UMAP", 
                                  colour_by = "dataset",
                                  point_alpha = 0.2, 
                                  add_ticks = FALSE,
                                  point_size = 0.75)
dr_bcell_sample$layers[[1]]$aes_params$colour <- NULL
dr_bcell_sample$layers[[1]]$aes_params$shape <- 16
dr_bcell_sample$layers[[1]]$mapping$colour <- dr_bcell_sample$layers[[1]]$mapping$fill
dr_bcell_sample <- dr_bcell_sample + 
  guides(colour = guide_legend(title = "Sample", override.aes = list(alpha = 1,
                                                                     size = 2), ncol = 2),
         fill = FALSE,
         shape = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$dataset)

# B cell marker figures

kappa_lambda_markers <- c("IGKC", "IGLC2", "IGLC3")
exprs <- logcounts(sce_bcell)[cellassign.utils::get_ensembl_id(kappa_lambda_markers, sce_bcell),]
expr_limits <- c(min(exprs), max(exprs))

sce_bcell_tmp <- sce_bcell

if (!is.null(args$winsorized_expression_threshold)) {
  logcounts(sce_bcell_tmp) <- pmin(logcounts(sce_bcell_tmp), args$winsorized_expression_threshold)
  expr_limits[2] <- min(expr_limits[2], args$winsorized_expression_threshold)
}

kappa_lambda_plots <- lapply(kappa_lambda_markers, function(mgene) {
  p <- plotReducedDim(sce_bcell_tmp,
                      use_dimred = "UMAP",
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_bcell_tmp),
                      point_alpha = 0.5,
                      point_size = 0.75,
                      add_ticks = FALSE)
  p$layers[[1]]$aes_params$colour <- NULL
  p$layers[[1]]$aes_params$shape <- 16
  p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$fill
  
  p <- p + 
    guides(fill = FALSE,
           colour = FALSE) + 
    xlab("UMAP-1") + 
    ylab("UMAP-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_fill_gradientn(colours = gradient_colours, 
                         limits = expr_limits) + 
    scale_colour_gradientn(colours = gradient_colours, 
                         limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})
names(kappa_lambda_plots) <- kappa_lambda_markers

kappa_lambda_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                        colour_mapping = gradient_colours,
                                                        legend_title = "Log normalized counts",
                                                        type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
kappa_lambda_legend <- cellassign.utils::extract_legend(kappa_lambda_legend)

# scvis figure
scvis_plot <- plotReducedDim(sce_scvis_merged, 
                             use_dimred = "scvis", 
                             colour_by = "category_full",
                             point_alpha = 0.1, 
                             add_ticks = FALSE,
                             point_size = 0.75) 

scvis_plot <- scvis_plot + 
  guides(fill = guide_legend(title = "Celltype", override.aes = list(alpha = 1,
                                                                     size = 2))) + 
  xlab("scvis-1") + 
  ylab("scvis-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature()

# T cell dr figure
dr_tcell <- plotReducedDim(sce_tcell, 
                           use_dimred = "UMAP", 
                           colour_by = "celltype_full",
                           point_alpha = 0.4, 
                           add_ticks = FALSE,
                           point_size = 0.75)
dr_tcell$layers[[1]]$aes_params$colour <- NULL
dr_tcell$layers[[1]]$aes_params$shape <- 16
dr_tcell$layers[[1]]$mapping$colour <- dr_tcell$layers[[1]]$mapping$fill
dr_tcell <- dr_tcell + 
  guides(colour = FALSE,
         fill = FALSE,
         shape = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$celltype)

dr_tcell_sample <- plotReducedDim(sce_tcell, 
                                  use_dimred = "UMAP", 
                                  colour_by = "dataset",
                                  point_alpha = 0.4, 
                                  add_ticks = FALSE,
                                  point_size = 0.75)
dr_tcell_sample$layers[[1]]$aes_params$colour <- NULL
dr_tcell_sample$layers[[1]]$aes_params$shape <- 16
dr_tcell_sample$layers[[1]]$mapping$colour <- dr_tcell_sample$layers[[1]]$mapping$fill
dr_tcell_sample <- dr_tcell_sample + 
  guides(colour = FALSE,
         fill = FALSE,
         shape = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$dataset)

# T and B cell subtype proportions
tcell_proportions_raw <- compute_celltype_proportions(sce_tcell, celltype_col = "celltype_full", suffix = "")
bcell_proportions_raw <- compute_celltype_proportions(sce_bcell %>% 
                                                        scater::filter(celltype_full %in% bcell_labels), 
                                                      celltype_col = "celltype_full", 
                                                      suffix = "")

meta <- colData(sce_tcell) %>%
  data.frame(check.names = FALSE) %>%
  dplyr::select(dataset, patient, timepoint) %>% 
  unique

tcell_proportions <- tcell_proportions_raw %>% 
  dplyr::left_join(meta) %>%
  reshape2::melt(id.vars = c("dataset", "patient", "timepoint"), measure.vars = tcell_labels,
                 variable.name = "celltype_full", value.name = "proportion")

bcell_proportions <- bcell_proportions_raw %>%
  dplyr::left_join(meta) %>%
  reshape2::melt(id.vars = c("dataset", "patient", "timepoint"), measure.vars = bcell_labels,
                 variable.name = "celltype_full", value.name = "proportion") %>%
  dplyr::mutate(celltype_full = plyr::mapvalues(celltype_full,
                                                from = c("B cells (malignant)", "B cells"),
                                                to = c("malignant", "nonmalignant"))) %>%
  dplyr::mutate(celltype_full = factor(celltype_full, levels = c("nonmalignant", "malignant")))

legend_limits <- round(c(0, max(c(tcell_proportions$proportion, bcell_proportions$proportion))), 1)

plot_proportions <- function(props, heatmap_heat_colours, 
                             legend_limits = NULL,
                             write_proportions = TRUE,
                             n_sig_dec_digits = 2,
                             plot_type = "lineplot",
                             celltype_palette,
                             legend = TRUE) {
  if (plot_type == "heatmap") {
    # Deprecated
    p <- ggplot(props, aes(x=timepoint, y=celltype_full)) + 
      geom_tile(aes(fill=proportion)) + 
      theme_bw() + 
      theme_Publication() + 
      theme_nature() + 
      scale_x_discrete(expand=c(0,0)) + 
      scale_y_discrete(expand=c(0,0)) + 
      theme(
        axis.line = element_blank()
      ) +
      xlab("Timepoint") + 
      ylab("Celltype") + 
      scale_fill_gradientn(colours = heatmap_heat_colours, limits = legend_limits) +
      theme(legend.position = "right", 
            legend.direction = "vertical", 
            legend.key.height = unit(4, "lines"),
            legend.key.width = unit(1, "lines"))
    if (legend) {
      p <- p + 
        guides(fill = guide_colorbar(title = "Proportion"))
    } else {
      p <- p + 
        guides(fill = FALSE)
    }
    
    if (write_proportions) {
      p <- p + 
        geom_text(aes(label=sapply(proportion, function(x) format(x, digits = n_sig_dec_digits, nsmall = n_sig_dec_digits))))
    }
  } else if (plot_type == "lineplot") {
    p <- ggplot(props, aes(x=timepoint, y=proportion, colour=celltype_full)) + 
      geom_point() + 
      geom_line(aes(group=celltype_full)) +
      theme_bw() + 
      theme_Publication() + 
      theme_nature() + 
      xlab("Timepoint") + 
      ylab("Proportion") + 
      scale_colour_manual(values = celltype_palette) 
    
    if (legend) {
      p <- p +
        guides(colour = guide_legend(title = "Celltype"))
    } else {
      p <- p + 
        guides(colour = FALSE)
    }
    
    p <- p + 
      facet_wrap(~ patient, nrow = 1) + 
      stripped_theme(strip_face = "bold")
    
  }
  
  return(p)
}

tcell_proportion_plot <- plot_proportions(tcell_proportions, 
                                          plot_type = "lineplot", 
                                          celltype_palette = categorical_palettes$celltype,
                                          legend = FALSE)

b_palette <- categorical_palettes$celltype[c("B cells", "B cells (malignant)")]
names(b_palette) <- names(b_palette) %>%
  plyr::mapvalues(from = c("B cells", "B cells (malignant)"),
                  to = c("nonmalignant", "malignant"))

bcell_proportion_plot <- plot_proportions(bcell_proportions, 
                                          plot_type = "lineplot", 
                                          celltype_palette = b_palette,
                                          legend = TRUE)


# CD8 T cell pathway/genes plot

de_cytotoxic_files <- Sys.glob(file.path(de_cytotoxic_dir, "*.rds"))
cytotoxic_de_results <- lapply(de_cytotoxic_files, function(f) {
  readRDS(f)
})
names(cytotoxic_de_results) <- tools::file_path_sans_ext(basename(de_cytotoxic_files))

# fl1018_cytotoxic_up <- cytotoxic_de_results$FL1018$pathway$up
# fl1018_cytotoxic_up@result <- fl1018_cytotoxic_up@result %>%
#   dplyr::mutate(old_description=Description,
#                 Description=)
# 
# emapplot(cytotoxic_de_results$FL1018$pathway$up, layout = 'kk')
# 
# plot_markers(cytotoxic_de_results$FL1018$gene, top_n = 10, top_direction = "equal",
#              highlight_genes = c("CD69", "IFNG"))

# CD8 T cell activation plot
azizi_signatures <- process_azizi_signatures(args$azizi_signatures)
cd8_activation_genes <- azizi_signatures$`CD8 T Cell Activation`
cd8_activation_genes <- cd8_activation_genes[!is.na(get_ensembl_id(cd8_activation_genes, sce_tcell))]
cd8_ensembl_ids <- get_ensembl_id(cd8_activation_genes, sce_tcell)

sce_cd8_activation <- sce_tcell[cd8_ensembl_ids,] %>%
  scater::filter(celltype_full %in% c("Cytotoxic T cells"))

selected_symbols <- c("GZMA", "PRF1", "CD69", "IFNG")
selected_ids <- get_ensembl_id(selected_symbols, sce_cd8_activation)

expr_table <- cbind(t(logcounts(sce_cd8_activation[selected_ids,])) %>% as.matrix %>% as.data.frame, 
                    colData(sce_cd8_activation))

expr_table_melted <- expr_table %>% 
  as.data.frame %>%
  reshape2::melt(id.vars = c("Sample", "Barcode", "dataset", "patient", "timepoint"),
                 measure.vars = selected_ids,
                 variable.name = "ID",
                 value.name = "expression") %>%
  dplyr::mutate(Symbol=df_as_map(rowData(sce_cd8_activation) %>% as.data.frame,
                                 ID,
                                 "ID",
                                 "Symbol")) %>%
  dplyr::mutate(Symbol=factor(Symbol, levels = selected_symbols))


cytotoxic_pvals <- compute_pvals_subsets(expr_table_melted, 
                                         facet_vars = c("patient", "Symbol"), 
                                         formula = expression ~ timepoint,
                                         corfun = wilcox.test,
                                         output = "p.value")
  
cytotoxic_boxplots <- ggplot(expr_table_melted, aes(x=timepoint, y=expression)) + 
  stat_boxplot(geom = 'errorbar', width = 0.3) + 
  geom_boxplot(outlier.size = -1, width = 0.5) + 
  geom_jitter(aes(colour=timepoint), position=position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  stripped_theme(strip_face = "bold") + 
  facet_grid(patient ~ Symbol) + 
  guides(colour = FALSE) + 
  xlab("Timepoint") + 
  ylab("Log normalized counts") + 
  geom_text(data = cytotoxic_pvals, aes(label=p.adj.text), 
            parse = TRUE, size = 0.35*8, x=Inf, y=Inf,
            hjust = 1.1, vjust = 1.3)



# Create a celltype-timepoint legend
tmp <- plotReducedDim(sce %>% scater::filter(celltype_full != "other"), use_dimred = args$dimreduce_type, colour_by = "celltype_full",
                      shape_by = "timepoint", point_alpha = 1) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$celltype) + 
  guides(colour=guide_legend(title = "Celltype"),
         shape=guide_legend(title = "Timepoint"))

celltype_timepoint_legend <- extract_legend(tmp)

# Create a T cell-specific celltype legend
tcell_legend <- cellassign.utils::ggsimplelegend(vars = names(categorical_palettes$celltype[tcell_labels]),
                                                 colour_mapping = categorical_palettes$celltype[tcell_labels],
                                                 legend_title = "Celltype",
                                                 legend_rows = 1)
tcell_legend <- cellassign.utils::extract_legend(tcell_legend)

sample_legend <- cellassign.utils::ggsimplelegend(vars = names(categorical_palettes$dataset[unique(sce_tcell$dataset)]),
                                                  colour_mapping = categorical_palettes$dataset[unique(sce_tcell$dataset)],
                                                  legend_title = "Sample",
                                                  legend_rows = 2)
sample_legend <- cellassign.utils::extract_legend(sample_legend)

# Build combined plot
kappa_lambda_plots_combined <- cowplot::plot_grid(plotlist = kappa_lambda_plots, nrow = 1)
kappa_lambda_plots_combined_withlegend <- cowplot::plot_grid(kappa_lambda_plots_combined,
                                                             kappa_lambda_legend,
                                                             nrow = 2,
                                                             rel_heights = c(0.8, 0.2))

top_row <- cowplot::plot_grid(dr_bcell_sample, 
                              dr_bcell,
                              kappa_lambda_plots_combined_withlegend, 
                              labels = c('a', '', 'b'),
                              ncol = 3,
                              rel_widths = c(0.25, 0.25, 0.5))

second_row <- cowplot::plot_grid(scvis_plot,
                                 bcell_proportion_plot,
                                 labels = c('c', 'd'),
                                 ncol = 2, 
                                 rel_widths = c(0.5, 0.5))


tcell_overall_combined <- cowplot::plot_grid(dr_tcell_sample,
                                             dr_tcell,
                                             tcell_proportion_plot,
                                             labels = c('e', '', 'f'),
                                             ncol = 3,
                                             rel_widths = c(0.25, 0.25, 0.5))
tcell_legend_row <- cowplot::plot_grid(sample_legend,
                                       tcell_legend,
                                       ncol = 2,
                                       rel_widths = c(0.25, 0.75))

third_row <- cowplot::plot_grid(tcell_overall_combined,
                                tcell_legend_row,
                                nrow = 2,
                                rel_heights = c(0.8, 0.2))

final_plot <- cowplot::plot_grid(top_row, 
                                 second_row, 
                                 third_row, 
                                 cytotoxic_boxplots,
                                 labels = c('', '', '', 'g'), 
                                 ncol = 1, 
                                 nrow = 4,
                                 rel_heights = c(0.3, 0.3, 0.3, 0.3))

# Plot final plot
pdf(args$outfname, width = 10, height = 13, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


