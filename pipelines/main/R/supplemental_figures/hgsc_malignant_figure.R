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
library(ggrastr)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create overview figure for HGSC")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--de_epithelial', metavar='DIR', type='character',
                    help="Directory to DE epithelial clusters results")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--hypoxia_winsorized_expression_threshold', type='double',
                    help="Hypoxia winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
de_epithelial_dir <- args$de_epithelial

categorical_palettes <- cat_palettes()
gradient_colours <- scrna_expression_gradient()

sce_epithelial <- sce %>%
  scater::filter(!is.na(epithelial_cluster))

# FGSEA plots

de_epithelial_fgsea_res <- readRDS(Sys.glob(file.path(de_epithelial_dir, "*.rds")))

fgsea_hla <- de_epithelial_fgsea_res$pathway %>%
  dplyr::filter(padj < 0.05,
                clust1 == "3",
                clust2 %in% c("1"))

fgsea_hypoxia <- de_epithelial_fgsea_res$pathway %>%
  dplyr::filter(padj < 0.05,
                clust1 == "2",
                clust2 %in% c("0", "4"))

nes_values <- c(fgsea_hla$NES,
                fgsea_hypoxia$NES)
nes_limits <- c(min(nes_values), max(nes_values))

fgsea_size_values <- c(fgsea_hla$size,
                       fgsea_hypoxia$size)
fgsea_size_limits <- c(min(fgsea_size_values), max(fgsea_size_values))

## HLA cluster-specific DE

fgsea_clust_31_plot <- ggplot(fgsea_hla %>%
                                dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                                              enriched_in=ifelse(NES < 0, 
                                                                 "1",
                                                                 "3")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=enriched_in)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,3)) + 
  scale_colour_manual(values = categorical_palettes$hgsc_epithelial_cluster) +
  guides(size = FALSE,
         colour = FALSE)

## HYPOXIA

fgsea_clust_20_plot <- ggplot(fgsea_hypoxia %>%
                                dplyr::filter(clust2 == "0") %>%
                                dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                                              enriched_in=ifelse(NES < 0, 
                                                                 "0",
                                                                 "2")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=enriched_in)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,3)) + 
  scale_colour_manual(values = categorical_palettes$hgsc_epithelial_cluster) +
  guides(size = FALSE,
         colour = FALSE)


fgsea_clust_24_plot <- ggplot(fgsea_hypoxia %>%
                                dplyr::filter(clust2 == "4") %>%
                                dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                                              enriched_in=ifelse(NES < 0, 
                                                                 "4",
                                                                 "2")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=enriched_in)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,3)) + 
  scale_colour_manual(values = categorical_palettes$hgsc_epithelial_cluster) +
  guides(size = FALSE,
         colour = FALSE)


# Volcano plot

highlight_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "HLA-DRA", "HLA-DRB1", "CD74")

gene_table <- de_epithelial_fgsea_res$gene %>%
  dplyr::filter(contrast == "3") %>%
  dplyr::select(logFDR.1, logFC.1, Symbol) %>%
  dplyr::mutate(FDR=exp(logFDR.1)) %>%
  dplyr::mutate(Significance=plyr::mapvalues(logFDR.1 < log(0.05),
                                             c(FALSE, TRUE), 
                                             c("P > 0.05", "P <= 0.05")))

xlims <- c(-1,1) * max(abs(gene_table$logFC.1))

hla_volcano_plot <- plot_markers_v2(gene_table,
                     top_n = 0,
                     top_direction = "equal",
                     highlight_genes = highlight_genes,
                     colour_col = "Significance", 
                     subtitle = "",
                     label_size = 2.5, 
                     point_alpha = 0.3,
                     force = 10,
                     min.segment.length = unit(0.1, 'lines')) + 
  scale_color_manual(values = categorical_palettes$significance) +
  guides(colour = FALSE) + 
  scale_x_continuous(limits = xlims)

# Hypoxia marker plots

hypoxia_marker_genes <- c("VEGFA", "CA9", "SLC2A1")

exprs <- logcounts(sce_epithelial)[cellassign.utils::get_ensembl_id(hypoxia_marker_genes, sce_epithelial),]
hypoxia_expr_limits <- c(min(exprs), max(exprs))

sce_tmp <- sce_epithelial

if (!is.null(args$hypoxia_winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$hypoxia_winsorized_expression_threshold)
  hypoxia_expr_limits[2] <- min(hypoxia_expr_limits[2], args$hypoxia_winsorized_expression_threshold)
}

hypoxia_marker_plots <- lapply(hypoxia_marker_genes, function(mgene) {
  p <- plotReducedDim(sce_tmp,
                      use_dimred = args$dimreduce_type,
                      colour_by = cellassign.utils::get_ensembl_id(mgene, sce_tmp),
                      point_alpha = 0.4,
                      point_size = 0.75)
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
                           limits = hypoxia_expr_limits) + 
    ggtitle(mgene)
  return(p)
})
names(hypoxia_marker_plots) <- hypoxia_marker_genes

# Legends

## Epithelial clusters

epithelial_cluster_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$hgsc_epithelial_cluster),
                                                              colour_mapping = categorical_palettes$hgsc_epithelial_cluster,
                                                              legend_title = "Cluster",
                                                              type = "discrete",
                                                              legend_rows = 1)
epithelial_cluster_legend <- cellassign.utils::extract_legend(epithelial_cluster_legend)

## Expression values

hypoxia_marker_legend <- cellassign.utils::ggsimplelegend(hypoxia_expr_limits,
                                                          colour_mapping = gradient_colours,
                                                          legend_title = "Log normalized counts",
                                                          type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
hypoxia_marker_legend <- cellassign.utils::extract_legend(hypoxia_marker_legend)

## fgsea legends

fgsea_colour_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$hgsc_epithelial_cluster),
                                                        colour_mapping = categorical_palettes$hgsc_epithelial_cluster,
                                                        legend_title = "Significantly upregulated in",
                                                        type = "discrete") 
fgsea_colour_legend <- cellassign.utils::extract_legend(fgsea_colour_legend)

fgsea_size_legend <- ggplot(fgsea_hla, aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,3)) + 
  guides(size = guide_legend(title = "Gene set size", title.position = "top", title.hjust = 0.5))
fgsea_size_legend <- cellassign.utils::extract_legend(fgsea_size_legend)

## Significance legend

significance_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$significance),
                                                        colour_mapping = categorical_palettes$significance,
                                                        legend_title = "Significance",
                                                        type = "discrete") 
significance_legend <- cellassign.utils::extract_legend(significance_legend)

# Assemble plot

hla_bottomrow <- cowplot::plot_grid(
  fgsea_clust_31_plot,
  hla_volcano_plot,
  nrow = 1,
  labels = c('a', 'b'))

fgsea_toprow <- cowplot::plot_grid(
  fgsea_clust_20_plot,
  fgsea_clust_24_plot,
  nrow = 1,
  labels = c('c', 'd')
)

legendrow <- cowplot::plot_grid(
  fgsea_colour_legend,
  fgsea_size_legend,
  significance_legend,
  nrow = 1
)

hypoxia_markerrow <- cowplot::plot_grid(
  plotlist = hypoxia_marker_plots,
  nrow = 1,
  labels = c('e', '', '', '')
)

final_plot <- cowplot::plot_grid(hla_bottomrow,
                                 legendrow,
                                 fgsea_toprow,
                                 hypoxia_markerrow,
                                 hypoxia_marker_legend,
                                 nrow = 5,
                                 rel_heights = c(0.8, 0.1, 1, 0.66, 0.1))



# Plot to output file
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")








