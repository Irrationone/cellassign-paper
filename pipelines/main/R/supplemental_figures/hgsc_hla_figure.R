# Overview figure showing HGSC cell types and HLA

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

parser <- ArgumentParser(description = "Create HGSC HLA figure")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--hla_winsorized_expression_threshold', type='double',
                    help="HLA winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

categorical_palettes <- cat_palettes()
gradient_colours <- scrna_expression_gradient()

sce_epithelial <- sce %>%
  scater::filter(!is.na(epithelial_cluster))


# HLA marker plots

hla_marker_genes <- c("HLA-A", "HLA-B", "HLA-C")

exprs <- logcounts(sce)[cellassign.utils::get_ensembl_id(hla_marker_genes, sce),]
hla_expr_limits <- c(min(exprs), max(exprs))

sce_tmp <- sce

if (!is.null(args$hla_winsorized_expression_threshold)) {
  logcounts(sce_tmp) <- pmin(logcounts(sce_tmp), args$hla_winsorized_expression_threshold)
  hla_expr_limits[2] <- min(hla_expr_limits[2], args$hla_winsorized_expression_threshold)
}

hla_marker_plots <- lapply(hla_marker_genes, function(mgene) {
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
                           limits = hla_expr_limits) + 
    ggtitle(mgene)
  return(p)
})
names(hla_marker_plots) <- hla_marker_genes

# HLA boxplots

hla_expr_cols <- exprs %>% as.matrix %>% t %>% as.data.frame
colnames(hla_expr_cols) <- hla_marker_genes

hla_expression_df <- colData(sce) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  cbind(hla_expr_cols) %>%
  reshape2::melt(measure.vars = hla_marker_genes, variable.name = "hla_gene", value.name = "expression")

hla_expression_df_filtered <- hla_expression_df %>%
  dplyr::filter(celltype != "other") %>%
  dplyr::mutate(group = ifelse(celltype == "Epithelial cells", ifelse(epithelial_cluster == "1", 
                                                                      "Epithelial (1)",
                                                                      "Epithelial (other)"), as.character(celltype)))
hla_mean_exprs <- hla_expression_df_filtered %>%
  dplyr::group_by(group) %>% 
  dplyr::summarise(mean_expr=mean(expression))

group_order <- (hla_mean_exprs %>% dplyr::arrange(-mean_expr))$group
hla_expression_df_filtered <- hla_expression_df_filtered %>%
  dplyr::mutate(group = factor(group, levels = group_order)) %>%
  dplyr::filter(!is.na(group))

# Almost everything is significant because of # of cells
# pvals <- setNames(plyr::ddply(hla_expression_df_filtered, plyr::.(hla_gene), function(x) {
#   df <- as.data.frame(x)
#   dunn.res <- dunn.test::dunn.test(df[,"expression"], df[,"group"], method = "bh", kw = TRUE, label = TRUE)
#   
#   return(data.frame(comparison=dunn.res$comparisons, p.adj=dunn.res$P.adjusted))
# }), c("hla_gene", "comparison", "p.value"))

hla_group_boxplots <- ggplot(hla_expression_df_filtered, aes(x=group, y=expression)) + 
  geom_boxplot(aes(fill=celltype), outlier.size = -1) + 
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.4, size = 1) +
  theme_bw() + 
  theme_Publication() +
  theme_nature() + 
  scale_fill_manual(values = categorical_palettes$hgsc_celltype) + 
  facet_wrap(~ hla_gene, ncol = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Group") + 
  ylab("Log normalized counts") + 
  stripped_theme(strip_face = "bold") + 
  guides(fill = FALSE)


# Legends


## Expression values

hla_marker_legend <- cellassign.utils::ggsimplelegend(hla_expr_limits,
                                                      colour_mapping = gradient_colours,
                                                      legend_title = "Log normalized counts",
                                                      type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
hla_marker_legend <- cellassign.utils::extract_legend(hla_marker_legend)


# Assemble plot

hla_markerrow <- cowplot::plot_grid(
  plotlist = hla_marker_plots,
  nrow = 1, 
  labels = c('a', '', '')
)


final_plot <- cowplot::plot_grid(hla_markerrow,
                                 hla_marker_legend,
                                 hla_group_boxplots,
                                 nrow = 3,
                                 labels = c('', '', 'b'),
                                 rel_heights = c(0.66, 0.08, 1))

pdf(args$outfname, width = 10, height = 9, useDingbats = FALSE)
plot(final_plot)
dev.off()


cat("Completed.\n")








