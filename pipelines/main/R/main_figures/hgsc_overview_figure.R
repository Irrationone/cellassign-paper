# Overview figure showing HGSC cell types

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)
library(scales)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create overview figure for HGSC")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_raw', metavar='FILE', type='character',
                    help="Path to raw SingleCellExperiment RDS")
parser$add_argument('--de_site', metavar='DIR', type='character',
                    help="Directory to DE site results")
parser$add_argument('--de_site_fgsea', metavar='DIR', type='character',
                    help="Directory to DE site fgsea results")
parser$add_argument('--de_epithelial', metavar='DIR', type='character',
                    help="Directory to DE epithelial clusters results")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--epithelial_winsorized_expression_threshold', type='double',
                    help="Epithelial winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce_raw_path <- args$sce_raw
sce <- readRDS(sce_path)
sce_raw <- readRDS(sce_raw_path)
de_site_dir <- args$de_site
de_site_fgsea_dir <- args$de_site_fgsea
de_epithelial_dir <- args$de_epithelial

categorical_palettes <- cat_palettes()

# Plot of timepoint
dr_site <- plotReducedDim(sce, use_dimred = args$dimreduce_type, 
                          colour_by = "dataset", point_alpha = 0.4, add_ticks = FALSE)
dr_site <- dr_site + 
  guides(colour = FALSE,
         fill = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$hgsc_dataset)
dr_site$layers[[1]]$aes_params$colour <- NULL
dr_site$layers[[1]]$aes_params$shape <- 16
dr_site$layers[[1]]$mapping$colour <- dr_site$layers[[1]]$mapping$fill

# Plot of celltype assignments
nonother_types <- sort(setdiff(unique(sce$celltype), "other"))
sce_celltype_remapped <- sce %>%
  scater::mutate(celltype=factor(plyr::mapvalues(celltype, from = c("other"),
                                                 to = c("Unassigned")),
                                 levels = c(nonother_types, "Unassigned")))
dr_celltype <- plotReducedDim(sce_celltype_remapped, 
                              use_dimred = args$dimreduce_type,
                              colour_by = "celltype",
                              point_alpha = 0.4, 
                              add_ticks = FALSE)
dr_celltype <- dr_celltype + 
  guides(colour = FALSE,
         fill = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$hgsc_celltype)
dr_celltype$layers[[1]]$aes_params$colour <- NULL
dr_celltype$layers[[1]]$aes_params$shape <- 16
dr_celltype$layers[[1]]$mapping$colour <- dr_celltype$layers[[1]]$mapping$fill

# Plots of marker gene expression
marker_genes <- c("EPCAM", "PTPRC", "MUM1L1", "COL1A1")
                  #"ACTA2","PECAM1", "MYH11", "MCAM")
exprs <- logcounts(sce)[cellassign.utils::get_ensembl_id(marker_genes, sce),]
expr_limits <- c(min(exprs), max(exprs))

gradient_colours <- scrna_expression_gradient()

sce_tmp <- sce

if (!is.null(args$winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$winsorized_expression_threshold)
  expr_limits[2] <- min(expr_limits[2], args$winsorized_expression_threshold)
}

marker_plots <- lapply(marker_genes, function(mgene) {
  p <- plotReducedDim(sce_tmp,
                      use_dimred = args$dimreduce_type,
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.2,
                      point_size = 0.5,
                      add_ticks = FALSE)
  p$layers[[1]]$aes_params$colour <- NULL
  p$layers[[1]]$aes_params$shape <- 16
  p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$fill
  
  p <- p + 
    guides(fill = FALSE,
           colour = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})

# Celltype proportion plot

nonother_types <- sort(setdiff(unique(sce$celltype), "other"))
coldat <- colData(sce_celltype_remapped) %>%
  as.data.frame

celltype_proportion_plot <- ggplot(coldat, aes(x=dataset)) + 
  geom_bar(aes(fill=celltype), position = "fill") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() +
  scale_y_continuous(labels = percent_format(), expand = c(0,0)) + 
  scale_fill_manual(values = categorical_palettes$hgsc_celltype) +
  xlab("Sample") + 
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = FALSE)


# T1 vs. T2 epithelial

de_site_fgsea_res <- readRDS(Sys.glob(file.path(de_site_fgsea_dir, "epithelial", "*.rds")))

epithelial_site_fgsea_filtered <- de_site_fgsea_res$pathway %>%
  dplyr::filter(padj < 0.05)

nes_values <- c(epithelial_site_fgsea_filtered$NES)
nes_limits <- c(min(nes_values), max(nes_values))

fgsea_size_values <- c(epithelial_site_fgsea_filtered$size)
fgsea_size_limits <- c(min(fgsea_size_values), max(fgsea_size_values))

fgsea_site_pathway_plot <- ggplot(epithelial_site_fgsea_filtered %>%
                                    dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                                                  NES=-1*NES,
                                                  enriched_in=ifelse(NES < 0, 
                                                                      "Right Ovary",
                                                                      "Left Ovary")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=enriched_in)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,4)) + 
  scale_colour_manual(values = categorical_palettes$hgsc_dataset) +
  guides(size = FALSE,
         colour = FALSE)

# Unsupervised clustering plot for epithelial clusters

sce_epithelial <- sce %>%
  scater::filter(!is.na(epithelial_cluster))

dr_epithelial_site <- plotReducedDim(sce_epithelial, use_dimred = args$dimreduce_type, 
                                     colour_by = "dataset", point_alpha = 0.5, add_ticks = FALSE)
dr_epithelial_site <- dr_epithelial_site + 
  guides(colour = FALSE,
         fill = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$hgsc_dataset)
dr_epithelial_site$layers[[1]]$aes_params$colour <- NULL
dr_epithelial_site$layers[[1]]$aes_params$shape <- 16
dr_epithelial_site$layers[[1]]$mapping$colour <- dr_epithelial_site$layers[[1]]$mapping$fill

# Plot of cluster assignments
dr_epithelial_cluster <- plotReducedDim(sce_epithelial,
                                        use_dimred = args$dimreduce_type,
                                        colour_by = "epithelial_cluster",
                                        point_alpha = 0.5, 
                                        add_ticks = FALSE)
dr_epithelial_cluster <- dr_epithelial_cluster + 
  guides(colour = FALSE,
         fill = FALSE) + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$hgsc_epithelial_cluster)
dr_epithelial_cluster$layers[[1]]$aes_params$colour <- NULL
dr_epithelial_cluster$layers[[1]]$aes_params$shape <- 16
dr_epithelial_cluster$layers[[1]]$mapping$colour <- dr_epithelial_cluster$layers[[1]]$mapping$fill

# Epithelial marker plots

epithelial_marker_genes <- c("CDH2", "THY1", "HLA-A", "HLA-B")

exprs <- logcounts(sce_epithelial)[cellassign.utils::get_ensembl_id(epithelial_marker_genes, sce_epithelial),]
epithelial_expr_limits <- c(min(exprs), max(exprs))

sce_tmp <- sce_epithelial

if (!is.null(args$epithelial_winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$epithelial_winsorized_expression_threshold)
  epithelial_expr_limits[2] <- min(epithelial_expr_limits[2], args$epithelial_winsorized_expression_threshold)
}

epithelial_marker_plots <- lapply(epithelial_marker_genes, function(mgene) {
  p <- plotReducedDim(sce_tmp,
                      use_dimred = args$dimreduce_type,
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.4,
                      point_size = 0.75,
                      add_ticks = FALSE)
  p$layers[[1]]$aes_params$colour <- NULL
  p$layers[[1]]$aes_params$shape <- 16
  p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$fill
  
  p <- p + 
    guides(fill = FALSE,
           colour = FALSE) + 
    xlab(paste0(args$dimreduce_type, "-1")) + 
    ylab(paste0(args$dimreduce_type, "-2")) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = epithelial_expr_limits) + 
    ggtitle(mgene)
  return(p)
})
names(epithelial_marker_plots) <- epithelial_marker_genes

# Legends

## Celltypes

celltype_legend <- cellassign.utils::ggsimplelegend(levels(sce_celltype_remapped$celltype),
                                                    colour_mapping = categorical_palettes$hgsc_celltype[levels(sce_celltype_remapped$celltype)],
                                                    legend_title = "Celltype",
                                                    type = "discrete",
                                                    legend_rows = 3)
celltype_legend <- cellassign.utils::extract_legend(celltype_legend)

## Sites

site_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$hgsc_dataset),
                                                colour_mapping = categorical_palettes$hgsc_dataset,
                                                legend_title = "Sample",
                                                type = "discrete",
                                                legend_rows = 1)
site_legend <- cellassign.utils::extract_legend(site_legend)

## Epithelial clusters

epithelial_cluster_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$hgsc_epithelial_cluster),
                                                              colour_mapping = categorical_palettes$hgsc_epithelial_cluster,
                                                              legend_title = "Cluster",
                                                              type = "discrete",
                                                              legend_rows = 1)
epithelial_cluster_legend <- cellassign.utils::extract_legend(epithelial_cluster_legend)

## Expression values

marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Log normalized counts",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)

epithelial_marker_legend <- cellassign.utils::ggsimplelegend(epithelial_expr_limits,
                                                             colour_mapping = gradient_colours,
                                                             legend_title = "Log normalized counts",
                                                             type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
epithelial_marker_legend <- cellassign.utils::extract_legend(epithelial_marker_legend)

## fgsea legends

fgsea_colour_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$hgsc_dataset),
                                                        colour_mapping = categorical_palettes$hgsc_dataset,
                                                        legend_title = "Significantly upregulated in",
                                                        type = "discrete") 
fgsea_colour_legend <- cellassign.utils::extract_legend(fgsea_colour_legend)

fgsea_size_legend <- ggplot(epithelial_site_fgsea_filtered, aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1, 4)) + 
  guides(size = guide_legend(title = "Gene set size", title.position = "top", title.hjust = 0.5))
fgsea_size_legend <- cellassign.utils::extract_legend(fgsea_size_legend)

# Assemble plot

## DR plots
dr_plots <- cowplot::plot_grid(dr_site, 
                               dr_celltype, 
                               celltype_proportion_plot, 
                               ncol = 3, 
                               nrow = 1, 
                               labels = c('a', 'b', 'c'),
                               rel_widths = c(0.4, 0.4, 0.2))
dr_plots_legend <- cowplot::plot_grid(site_legend, celltype_legend,
                                      ncol = 2,
                                      rel_widths = c(0.4, 0.6))
marker_gene_plots <- cowplot::plot_grid(plotlist = marker_plots, ncol = 4, labels = c('d', '', '', ''))


epithelial_plots <- cowplot::plot_grid(fgsea_site_pathway_plot,
                                       dr_epithelial_site,
                                       dr_epithelial_cluster,
                                       ncol = 3,
                                       labels = c('e', 'f', ''),
                                       rel_widths = c(0.5, 0.25, 0.25))
epithelial_plots_legend <- cowplot::plot_grid(fgsea_size_legend,
                                              fgsea_colour_legend,
                                              site_legend,
                                              epithelial_cluster_legend,
                                              ncol = 4)

epithelial_marker_plots_combined <- cowplot::plot_grid(plotlist = epithelial_marker_plots,
                                                       ncol = 4,
                                                       labels = c('g', '', 'h', ''))

final_plot <- cowplot::plot_grid(dr_plots, 
                                 dr_plots_legend,
                                 marker_gene_plots, 
                                 marker_legend,
                                 epithelial_plots,
                                 epithelial_plots_legend,
                                 epithelial_marker_plots_combined,
                                 epithelial_marker_legend,
                                 labels = c('', '', ''), 
                                 ncol = 1, 
                                 nrow = 8,
                                 rel_heights = c(0.7, 0.1, 0.5, 0.1, 0.7, 0.1, 0.5, 0.1))

# Plot to output file
pdf(args$outfname, width = 10, height = 14, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")



