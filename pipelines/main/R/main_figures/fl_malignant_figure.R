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
library(org.Hs.eg.db)
library(ReactomePA)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create nonmalignant figure for FL")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_bcell', metavar='FILE', type='character',
                    help="Path to B cell-filtered SCE")
parser$add_argument('--bcell_labels', type='character', nargs='+',
                    help="Cell type labels of B cells")
parser$add_argument('--de_timepoint_dir', type='character',
                    help="DE results for timepoint comparisons")
parser$add_argument('--de_timepoint_fgsea_dir', type='character',
                    help="DE results for timepoint comparisons (using FGSEA and hallmark gene set)")
parser$add_argument('--de_malignant_timepoint_dir', type='character',
                    help="DE results for malignant-timepoint interaction comparisons")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
sce_bcell <- readRDS(args$sce_bcell)
bcell_labels <- unlist(args$bcell_labels)
sce_bcell <- sce_bcell %>%
  scater::filter(celltype_full %in% bcell_labels)

de_timepoint_dir <- args$de_timepoint_dir
de_timepoint_fgsea_dir <- args$de_timepoint_fgsea_dir
de_malignant_timepoint_dir <- args$de_malignant_timepoint_dir

de_timepoint_files <- Sys.glob(file.path(de_timepoint_dir, "malignant", "*"))
de_timepoint_fgsea_files <- Sys.glob(file.path(de_timepoint_fgsea_dir, "malignant", "*"))
de_malignant_timepoint_files <- Sys.glob(file.path(de_malignant_timepoint_dir, "*", "*"))


categorical_palettes <- cat_palettes()
heatmap_heat_colours <- heat_colour_gradient()
gradient_colours <- scrna_expression_gradient()
pval_colours <- pval_colour_gradient()

# Proliferation markers & cell cycle plots

patients <- sce_bcell$patient %>% unique
proliferation_markers <- c("MKI67", "TOP2A")
exprs <- logcounts(sce_bcell)[cellassign.utils::get_ensembl_id(proliferation_markers, sce_bcell),]
expr_limits <- c(min(exprs), max(exprs))

sce_bcell_tmp <- sce_bcell

if (!is.null(args$winsorized_expression_threshold)) {
  logcounts(sce_bcell_tmp) <- pmin(logcounts(sce_bcell_tmp), args$winsorized_expression_threshold)
  expr_limits[2] <- min(expr_limits[2], args$winsorized_expression_threshold)
}

proliferation_plots <- lapply(patients, function(pat) {
  res <- lapply(proliferation_markers, function(mgene) {
    p <- plotReducedDim(sce_bcell_tmp %>% scater::filter(malignant_status_manual == "malignant",
                                                         patient == pat),
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
  names(res) <- proliferation_markers
  return(res)
})
names(proliferation_plots) <- patients


# B cell dr figures
bcell_timepoint_plots <- lapply(patients, function(pat) {
  bcell_timepoint <- plotReducedDim(sce_bcell %>% scater::filter(malignant_status_manual == "malignant",
                                                                 patient == pat), 
                                    use_dimred = "UMAP", 
                                    colour_by = "timepoint",
                                    point_alpha = 0.2, 
                                    add_ticks = FALSE,
                                    point_size = 0.75)
  bcell_timepoint$layers[[1]]$aes_params$colour <- NULL
  bcell_timepoint$layers[[1]]$aes_params$shape <- 16
  bcell_timepoint$layers[[1]]$mapping$colour <- bcell_timepoint$layers[[1]]$mapping$fill
  
  bcell_timepoint <- bcell_timepoint + 
    guides(colour = FALSE,
           fill = FALSE,
           shape = FALSE) + 
    xlab("UMAP-1") + 
    ylab("UMAP-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() + 
    scale_colour_manual(values = categorical_palettes$timepoint)
  return(bcell_timepoint)
})
names(bcell_timepoint_plots) <- patients

# Cell cycle pairplots

cycling_stats <- with(colData(sce), table(patient, timepoint, celltype_full, Cell_Cycle)) %>% 
  as.data.frame %>%
  dplyr::rename(count=Freq) %>%
  dplyr::group_by(patient, timepoint, celltype_full) %>% 
  dplyr::mutate(total_count=sum(count)) %>%
  dplyr::ungroup()

cycling_summary <- cycling_stats %>% 
  dplyr::filter(Cell_Cycle %in% c("S", "G2M")) %>% 
  dplyr::group_by(patient, timepoint, celltype_full, total_count) %>%
  dplyr::summarise(cycling_count=sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cycling_prop=cycling_count/total_count) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(total_prop=total_count/sum(total_count)) %>%
  dplyr::ungroup()

min_cell_count <- 15

prop_limits <- c(0.01, max(cycling_summary$total_prop))
cycling_plot_results <- lapply(patients, function(pat) {
  cycling_summary_patient <- cycling_summary %>% 
    dplyr::filter(patient == pat)
  
  removed_cell_types <- (cycling_summary_patient %>% 
    dplyr::filter(total_count <= min_cell_count))$celltype_full %>% 
    unique %>%
    union(c("other"))
  
  included_cell_types <- cycling_summary_patient$celltype_full %>% 
    unique %>%
    setdiff(removed_cell_types)
  
  cycling_plot <- ggplot(cycling_summary_patient %>%
                           dplyr::filter(!celltype_full %in% removed_cell_types),
                         aes(x=timepoint, y=cycling_prop, colour=celltype_full)) +
    geom_point(aes(size=total_prop)) +
    geom_line(aes(group=celltype_full)) + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    xlab("Timepoint") + 
    ylab("Proportion of cells in S/G2/M") + 
    scale_size_continuous(trans = "log", breaks = c(0.01, 0.1, 0.5), limits = prop_limits) + 
    scale_colour_manual(values = categorical_palettes$celltype) + 
    guides(colour= FALSE,
           size = FALSE)
  
  return(list(plot=cycling_plot, cell_types=included_cell_types))
})
names(cycling_plot_results) <- patients

# Pathway plots

de_timepoint_res <- lapply(de_timepoint_files, function(f) {
  readRDS(f)
})
names(de_timepoint_res) <- tools::file_path_sans_ext(basename(de_timepoint_files))

de_timepoint_fgsea_res <- lapply(de_timepoint_fgsea_files, function(f) {
  readRDS(f)
})
names(de_timepoint_fgsea_res) <- tools::file_path_sans_ext(basename(de_timepoint_fgsea_files))

de_malignant_timepoint_res <- lapply(de_malignant_timepoint_files, function(f) {
  readRDS(f)
})
names(de_malignant_timepoint_res) <- tools::file_path_sans_ext(basename(de_malignant_timepoint_files))

FL1018_de_malignant_down_table <- de_malignant_timepoint_res$FL1018$malignant$pathway$down
FL1018_de_malignant_timepoint_down_table <- de_timepoint_res$FL1018$pathway$down

max_pathways <- 10
padj_values <- c(FL1018_de_malignant_down_table@result[1:max_pathways,"p.adjust"],
                 FL1018_de_malignant_timepoint_down_table@result[1:max_pathways,"p.adjust"])
padj_limits <- log10(c(min(padj_values), max(padj_values)))
padj_limits <- 10^(c(floor(padj_limits[1]), ceiling(padj_limits[2])))

sizes <- c(FL1018_de_malignant_down_table@result[1:max_pathways,"Count"],
           FL1018_de_malignant_timepoint_down_table@result[1:max_pathways,"Count"])
size_limits <- c(min(sizes), max(sizes))

## Malignant results for the INTERACTION model
FL1018_malignant_pathway_down_plot <- emapplot(FL1018_de_malignant_down_table,
                                               showCategory = max_pathways,
                                               layout = 'kk') + 
  guides(colour = FALSE,
         size = FALSE) + 
  scale_colour_gradientn(colours = rev(pval_colours),
                         trans = "log10",
                         limits = padj_limits) + 
  scale_size_continuous(limits = size_limits)

FL1018_malignant_pathway_down_plot$layers[[1]]$aes_params$edge_alpha <- 0.3
FL1018_malignant_pathway_down_plot$layers[[1]]$aes_params$edge_width <- 1
FL1018_malignant_pathway_down_plot$layers[[3]]$aes_params$size <- 2.5

## Timepoint results for malignant cells between T1 and T2. NOT the interaction model
FL1018_malignant_T2_pathway_down_plot <- emapplot(FL1018_de_malignant_timepoint_down_table,
                                                  showCategory = max_pathways,
                                                  layout = 'kk') + 
  guides(colour = FALSE,
         size = FALSE)  + 
  scale_colour_gradientn(colours = rev(pval_colours),
                         trans = "log10",
                         limits = padj_limits)  + 
  scale_size_continuous(limits = size_limits)

FL1018_malignant_T2_pathway_down_plot$layers[[1]]$aes_params$edge_alpha <- 0.3
FL1018_malignant_T2_pathway_down_plot$layers[[1]]$aes_params$edge_width <- 1
FL1018_malignant_T2_pathway_down_plot$layers[[3]]$aes_params$size <- 2.5


# Pathway plots for malignant over time

# Always include cell cycle pathways
force_include_pathways <- c("HALLMARK_MITOTIC_SPINDLE",
                            "HALLMARK_E2F_TARGETS",
                            "HALLMARK_G2M_CHECKPOINT")

FL1018_fgsea_filtered <- de_timepoint_fgsea_res$FL1018$pathway %>% 
  dplyr::filter(padj < 0.05 | pathway %in% force_include_pathways)
FL2001_fgsea_filtered <- de_timepoint_fgsea_res$FL2001$pathway %>% 
  dplyr::filter(padj < 0.05 | pathway %in% force_include_pathways)

nes_values <- c(FL1018_fgsea_filtered$NES,
                FL2001_fgsea_filtered$NES)
nes_limits <- c(min(nes_values), max(nes_values))

fgsea_size_values <- c(FL1018_fgsea_filtered$size,
                       FL2001_fgsea_filtered$size)
fgsea_size_limits <- c(min(fgsea_size_values), max(fgsea_size_values))

fgsea_pathway_plot_FL1018 <- ggplot(FL1018_fgsea_filtered %>%
         dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                       significance=ifelse(padj <= 0.05, 
                                           "P <= 0.05",
                                           "P > 0.05")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=significance)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,4)) + 
  scale_colour_manual(values = categorical_palettes$significance) +
  guides(size = FALSE,
         colour = FALSE)

fgsea_pathway_plot_FL2001 <- ggplot(FL2001_fgsea_filtered %>%
                                      dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", ""),
                                                    significance=ifelse(padj <= 0.05, 
                                                                        "P <= 0.05",
                                                                        "P > 0.05")), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=significance)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_y_continuous(limits = nes_limits) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  scale_size_continuous(trans = "log10", limits = fgsea_size_limits,
                        range = c(1,4)) + 
  scale_colour_manual(values = categorical_palettes$significance) +
  guides(size = FALSE,
         colour = FALSE)


# HLA expression plots

hla_genes <- c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1")
hla_ids <- get_ensembl_id(hla_genes, sce_bcell)

hla_exprs <- logcounts(sce_bcell)[hla_ids,] %>%
  as.matrix %>%
  t
colnames(hla_exprs) <- hla_genes
hla_exprs_annotated <- colData(sce_bcell) %>%
  data.frame(check.names = FALSE) %>%
  cbind(hla_exprs) %>%
  reshape2::melt(
    id.vars = c("patient", "timepoint", "dataset", "malignant_status_manual"),
    measure.vars = hla_genes,
    variable.name = "Symbol",
    value.name = "logcounts"
  ) %>%
  dplyr::mutate(timepoint = factor(timepoint),
                Symbol = factor(Symbol),
                malignant_status_manual = factor(malignant_status_manual, levels = c("nonmalignant", "malignant")))

hla_exprs_annotated_summary <- hla_exprs_annotated %>%
  dplyr::group_by(patient, dataset, timepoint, malignant_status_manual, Symbol) %>%
  dplyr::summarise(logcount_mean=mean(logcounts)) %>%
  dplyr::ungroup()

hla_boxplots_faceted <- ggplot(hla_exprs_annotated %>%
                                 dplyr::mutate(xlab=paste(malignant_status_manual, timepoint, sep = ", ")),
                               aes(x=xlab, y = logcounts)) +
  geom_boxplot(aes(fill=xlab)) +
  theme_bw() + 
  theme_Publication() +
  theme_nature() + 
  stripped_theme(strip_face = "bold") +
  facet_grid(patient ~ Symbol, scales = "free") +
  xlab("B cell population") +
  ylab("Log normalized counts") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  scale_fill_manual(values = categorical_palettes$bcell_population) +
  guides(fill = FALSE)

# Create legends

## Timepoint and expression legends

timepoint_legend <- cellassign.utils::ggsimplelegend(vars = names(categorical_palettes$timepoint),
                                                     colour_mapping = categorical_palettes$timepoint,
                                                     legend_title = "Timepoint",
                                                     legend_rows = 1)
timepoint_legend <- cellassign.utils::extract_legend(timepoint_legend)

proliferation_expression_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                                    colour_mapping = gradient_colours,
                                                                    legend_title = "Expression",
                                                                    type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
proliferation_expression_legend <- cellassign.utils::extract_legend(proliferation_expression_legend)

## Cell type and size legends for cell cycle plots

cycling_cell_types <- do.call(union, unname(lapply(cycling_plot_results, function(x) x$cell_types)))

cycling_celltype_legend <- cellassign.utils::ggsimplelegend(vars = names(categorical_palettes$celltype[cycling_cell_types]),
                                                            colour_mapping = categorical_palettes$celltype[cycling_cell_types],
                                                            legend_title = "Celltype",
                                                            legend_rows = 2)
cycling_celltype_legend <- cellassign.utils::extract_legend(cycling_celltype_legend)

message("Creating cycling size legend ...")

cycling_size_legend <- ggplot(cycling_summary, aes(x=timepoint, y=cycling_prop)) +
  geom_point(aes(size=total_prop)) +
  theme_bw() + 
  theme_Publication() + 
  theme_nature() +
  guides(size = guide_legend(title = "Total proportion",
                             title.position = "top",
                             title.hjust = 0.5)) + 
  scale_size_continuous(trans = "log", breaks = c(0.01, 0.1, 0.3), limits = prop_limits) + 
  theme(legend.title = element_text(size = 7))

cycling_size_legend <- cellassign.utils::extract_legend(cycling_size_legend)

message("Completed cycling size plot.")

# HLA legend

# hla_gene_legend <- cellassign.utils::ggsimplelegend(vars = names(categorical_palettes$hla_genes[hla_genes]),
#                                                     colour_mapping = categorical_palettes$hla_genes[hla_genes],
#                                                     legend_title = "Gene",
#                                                     legend_rows = 1)
# hla_gene_legend <- cellassign.utils::extract_legend(hla_gene_legend)

# Network legends

network_pval_legend <- cellassign.utils::ggsimplelegend(padj_limits,
                                                        colour_mapping = rev(pval_colours),
                                                        legend_title = "Adjusted p value",
                                                        type = "continuous") + 
  theme(legend.key.width = unit(2, "lines")) + 
  scale_fill_gradientn(colours = rev(pval_colours), limits = padj_limits, trans = "log10")
network_pval_legend <- cellassign.utils::extract_legend(network_pval_legend)

message("Creating network size legend ...")

network_size_legend <- emapplot(FL1018_de_malignant_timepoint_down_table,
                                showCategory = max_pathways,
                                layout = 'kk') + 
  guides(colour = FALSE,
         size = guide_legend(title = "Genes", 
                             title.position = "top", 
                             title.hjust = 0.5))  + 
  scale_colour_gradientn(colours = rev(pval_colours),
                         trans = "log10",
                         limits = padj_limits)  + 
  scale_size_continuous(limits = size_limits, 
                        breaks = c(3, 6, 12)) + 
  theme_nature(fontsize = 7) + 
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.background = element_rect(fill = "transparent", 
                                         colour = NA), 
        legend.justification = "center") 
network_size_legend <- cellassign.utils::extract_legend(network_size_legend)

# fgsea size legend

fgsea_size_legend <- ggplot(de_timepoint_fgsea_res$FL2001$pathway %>% 
                              dplyr::filter(padj < 0.05) %>%
                              dplyr::mutate(pathway=str_replace_all(pathway, "^HALLMARK_", "")), aes(reorder(pathway, NES), NES)) +
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
  guides(size = guide_legend(title = "Gene set size"))
fgsea_size_legend <- cellassign.utils::extract_legend(fgsea_size_legend)


fgsea_colour_legend <- cellassign.utils::ggsimplelegend(names(categorical_palettes$significance),
                                                        colour_mapping = categorical_palettes$significance,
                                                        legend_title = "Significance",
                                                        type = "discrete") 
fgsea_colour_legend <- cellassign.utils::extract_legend(fgsea_colour_legend)

# Build combined plot

fgsea_pathway_plot_FL1018 <- fgsea_pathway_plot_FL1018 %>%
  ggplot_build() %>%
  ggplot_gtable()

fgsea_pathway_plot_FL2001 <- fgsea_pathway_plot_FL2001 %>%
  ggplot_build() %>%
  ggplot_gtable()

fgsea_pathway_plot_FL1018$layout$clip[str_detect(fgsea_pathway_plot_FL1018$layout$name, "panel")] <- "off"
fgsea_pathway_plot_FL2001$layout$clip[str_detect(fgsea_pathway_plot_FL2001$layout$name, "panel")] <- "off"


fgsea_row <- cowplot::plot_grid(fgsea_pathway_plot_FL1018,
                                fgsea_pathway_plot_FL2001,
                                labels = c('a', 'b'),
                                ncol = 2)

fgsea_legend_row <- cowplot::plot_grid(fgsea_size_legend, 
                                       fgsea_colour_legend,
                                       labels = c('', ''),
                                       ncol = 2)

cycling_plot_results$FL1018$plot <- cycling_plot_results$FL1018$plot %>%
  ggplot_build() %>%
  ggplot_gtable()
cycling_plot_results$FL1018$plot$layout$clip[str_detect(cycling_plot_results$FL1018$plot$layout$name, "panel")] <- "off"

cycling_plot_results$FL2001$plot <- cycling_plot_results$FL2001$plot %>%
  ggplot_build() %>%
  ggplot_gtable()
cycling_plot_results$FL2001$plot$layout$clip[str_detect(cycling_plot_results$FL2001$plot$layout$name, "panel")] <- "off"


top_row <- cowplot::plot_grid(bcell_timepoint_plots$FL1018,
                              proliferation_plots$FL1018$MKI67,
                              proliferation_plots$FL1018$TOP2A,
                              cycling_plot_results$FL1018$plot,
                              labels = c('c', '', '', 'd'),
                              ncol = 4,
                              rel_widths = rep(0.25, 4))

second_row <- cowplot::plot_grid(bcell_timepoint_plots$FL2001,
                                 proliferation_plots$FL2001$MKI67,
                                 proliferation_plots$FL2001$TOP2A,
                                 cycling_plot_results$FL2001$plot,
                                 labels = c('e', '', '', 'f'),
                                 ncol = 4,
                                 rel_widths = rep(0.25, 4))

top_second_legend <- cowplot::plot_grid(timepoint_legend,
                                        proliferation_expression_legend,
                                        cycling_celltype_legend,
                                        cycling_size_legend, 
                                        ncol = 4,
                                        rel_widths = c(0.2, 0.35, 0.25, 0.20))

# third_row <- cowplot::plot_grid(FL1018_malignant_pathway_down_plot,
#                                 FL1018_malignant_T2_pathway_down_plot,
#                                 labels = c('e', 'f'),
#                                 ncol = 2,
#                                 rel_widths = c(0.5, 0.5))
# 
# third_row_legend <- cowplot::plot_grid(network_pval_legend,
#                                        network_size_legend,
#                                        ncol = 2,
#                                        rel_widths = c(0.5, 0.5))

final_plot <- cowplot::plot_grid(fgsea_row,
                                 fgsea_legend_row,
                                 top_row, 
                                 top_second_legend,
                                 second_row, 
                                 #third_row, 
                                 #third_row_legend,
                                 hla_boxplots_faceted,
                                 labels = c('', '', '', '', '', 'g'), 
                                 ncol = 1, 
                                 nrow = 6,
                                 rel_heights = c(0.4, 0.05, 0.25, 0.05, 0.25, 0.4))

# Plot final plot
pdf(args$outfname, width = 10, height = 13, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


