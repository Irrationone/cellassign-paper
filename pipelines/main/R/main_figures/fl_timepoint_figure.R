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

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create overview figure for FL")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_tcell', metavar='FILE', type='character',
                    help="Path to T cell-filtered SCE")
parser$add_argument('--sce_bcell', metavar='FILE', type='character',
                    help="Path to B cell-filtered SCE")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--tcell_labels', type='character', nargs='+',
                    help="Cell type labels of T cells")
parser$add_argument('--bcell_labels', type='character', nargs='+',
                    help="Cell type labels of B cells")
parser$add_argument('--azizi_signatures', metavar = 'FILE', type='character',
                    help="XLS file of Azizi et al. Table S4")
parser$add_argument('--de_malignant', type='character',
                    help="DE results for malignant B cells")
parser$add_argument('--de_nonmalignant', type='character',
                    help="DE results for nonmalignant B cells")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
sce_tcell <- readRDS(args$sce_tcell)
sce_bcell <- readRDS(args$sce_bcell)
de_malignant <- readRDS(args$de_malignant)
de_nonmalignant <- readRDS(args$de_nonmalignant)

tcell_labels <- unlist(args$tcell_labels)
bcell_labels <- unlist(args$bcell_labels)

categorical_palettes <- cat_palettes()
heatmap_heat_colours <- heat_colour_gradient()

# Process T and B cell subset SCEs
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
  

# T cell dr figure
dr_tcell <- plotReducedDim(sce_tcell, 
                           use_dimred = "UMAP", 
                           colour_by = "celltype_full",
                           shape_by = "timepoint",
                           point_alpha = 0.5, 
                           point_size = 2)
dr_tcell <- dr_tcell + 
  geom_rug(alpha = 0.1, colour = "gray20") +
  guides(colour = FALSE,
         shape = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_fill_manual(values = categorical_palettes$celltype)

# T and B cell subtype proportions
tcell_proportions_raw <- compute_celltype_proportions(sce_tcell, celltype_col = "celltype_full", suffix = "")
bcell_proportions_raw <- compute_celltype_proportions(sce_bcell %>% 
                                                        scater::filter(celltype_full %in% bcell_labels), 
                                                      celltype_col = "celltype_full", 
                                                      suffix = "")

meta <- colData(sce_tcell) %>%
  data.frame(check.names = FALSE) %>%
  dplyr::select(dataset, timepoint) %>% 
  unique

tcell_proportions <- tcell_proportions_raw %>% 
  dplyr::left_join(meta) %>%
  reshape2::melt(id.vars = c("dataset", "timepoint"), measure.vars = tcell_labels,
                 variable.name = "celltype_full", value.name = "proportion")

bcell_proportions <- bcell_proportions_raw %>%
  dplyr::left_join(meta) %>%
  reshape2::melt(id.vars = c("dataset", "timepoint"), measure.vars = bcell_labels,
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
    
  }
  
  return(p)
}

tcell_proportion_plot <- plot_proportions(tcell_proportions, 
                                          plot_type = "lineplot", 
                                          celltype_palette = categorical_palettes$celltype,
                                          legend = FALSE) +
  ggtitle("T cells")

b_palette <- categorical_palettes$celltype[c("B cells", "B cells (malignant)")]
names(b_palette) <- names(b_palette) %>%
  plyr::mapvalues(from = c("B cells", "B cells (malignant)"),
                  to = c("nonmalignant", "malignant"))

bcell_proportion_plot <- plot_proportions(bcell_proportions, 
                                          plot_type = "lineplot", 
                                          celltype_palette = b_palette,
                                          legend = FALSE) +
  ggtitle("B cells")


# CD8 T cell activation plot
azizi_signatures <- process_azizi_signatures(args$azizi_signatures)
cd8_activation_genes <- azizi_signatures$`CD8 T Cell Activation`
cd8_activation_genes <- cd8_activation_genes[!is.na(get_ensembl_id(cd8_activation_genes, sce_tcell))]
cd8_ensembl_ids <- get_ensembl_id(cd8_activation_genes, sce_tcell)

sce_cd8_activation <- sce_tcell[cd8_ensembl_ids,]
sce_t1_cytotoxic <- sce_cd8_activation %>%
  scater::filter(timepoint == "T1" & celltype_full %in% c("Cytotoxic T cells",
                                                          "Cytotoxic T cells (activated)"))
sce_t2_cytotoxic <- sce_cd8_activation %>%
  scater::filter(timepoint == "T2" & celltype_full %in% c("Cytotoxic T cells",
                                                          "Cytotoxic T cells (activated)"))
t1_cd8_activation <- scater::calcAverage(sce_t1_cytotoxic,
                                         exprs_values = "logcounts",
                                         use_size_factors = rep(1, ncol(sce_cd8_activation)))
t2_cd8_activation <- scater::calcAverage(sce_t2_cytotoxic,
                                         exprs_values = "logcounts",
                                         use_size_factors = rep(1, ncol(sce_cd8_activation)))

cd8_activation_scores <- rbind(t1_cd8_activation, t2_cd8_activation)
rownames(cd8_activation_scores) <- c("T1", "T2")
colnames(cd8_activation_scores) <- cd8_activation_genes

cd8_activation_scores_melted <- cd8_activation_scores %>%
  reshape2::melt() %>% 
  dplyr::rename(Timepoint=Var1,
                Symbol=Var2,
                mean_logcounts=value)


cd8_activation_scores_cells <- rbind(data.frame(timepoint='T1', mean_logcounts_genes=colMeans(logcounts(sce_t1_cytotoxic))),
                                     data.frame(timepoint='T2', mean_logcounts_genes=colMeans(logcounts(sce_t2_cytotoxic))))

pvals <- compute_pvals_subsets(cd8_activation_scores_cells %>%
                                 dplyr::mutate(temp=1), 
                               facet_vars = c("temp"), 
                               formula = mean_logcounts_genes ~ timepoint,
                               corfun = wilcox.test)

## CD8 activation boxplot
cd8_activation_box <- ggplot(cd8_activation_scores_cells, aes(x=timepoint, y=mean_logcounts_genes)) +
  geom_boxplot(width = 0.5, outlier.size = -1) + 
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  coord_flip() + 
  xlab("Timepoint") + 
  ylab("CD8 Activation") + 
  scale_x_discrete(limits = rev(levels(cd8_activation_scores_cells$timepoint)), expand = c(0.25,0.25)) + 
  geom_text(data = pvals, aes(x=Inf, y=Inf, label = p.adj.text), hjust = 1, vjust = 1, parse = TRUE, size = (0.35/1) * 8)

cd8_activation_box <- annotate_boxplot_p(cd8_activation_box,
                                         cd8_activation_scores_cells %>% dplyr::mutate(temp=1),
                                         facet_vars = c("temp"),
                                         tip_vars = c("temp"),
                                         x_name = "timepoint",
                                         y_name = "mean_logcounts_genes",
                                         show_threshold = 0.05,
                                         pvals = pvals %>% dplyr::mutate(comparison = "T2 - T1"),
                                         angle = 90,
                                         seps = 70)

gene_hc <- hclust(dist(t(cd8_activation_scores)), method = "ward.D2")

## CD8 activation heatmap
cd8_activation_heatmap <- ggplot(cd8_activation_scores_melted, aes(y=Timepoint, x=Symbol)) +
  geom_tile(aes(fill=mean_logcounts)) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_x_discrete(expand=c(0,0), limits = levels(cd8_activation_scores_melted$Symbol)[gene_hc$order]) + 
  scale_y_discrete(expand=c(0,0), limits = rev(levels(cd8_activation_scores_melted$Timepoint))) + 
  theme(
    axis.line = element_blank()
  ) +
  xlab("Gene") + 
  ylab("") + 
  scale_fill_gradientn(colours = heatmap_heat_colours, limits = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  guides(fill = FALSE) 

cd8_activation_plot <- gridExtra::gtable_cbind(ggplotGrob(cd8_activation_box), 
                                               ggplotGrob(cd8_activation_heatmap))


# HLA downregulation in malignant B cells T2 vs. T1
hla_genes <- c("B2M", "HLA-A", "HLA-B", "HLA-C")#, "CD74", "HLA-DRA", "HLA-DRB1")
hla_ids <- get_ensembl_id(hla_genes, sce)
malignant_de_genes <- de_malignant$genes %>%
  data.frame(check.names = FALSE)
nonmalignant_de_genes <- de_nonmalignant$genes %>%
  data.frame(check.names = FALSE)

malignant_de_hla <- malignant_de_genes %>% 
  dplyr::filter(Symbol %in% hla_genes)
nonmalignant_de_hla <- nonmalignant_de_genes %>%
  dplyr::filter(Symbol %in% hla_genes)
hla_exprs <- logcounts(sce_bcell)[hla_ids,] %>%
  as.matrix %>%
  t
colnames(hla_exprs) <- hla_genes
hla_exprs_annotated <- colData(sce_bcell) %>%
  data.frame(check.names = FALSE) %>%
  cbind(hla_exprs) %>%
  reshape2::melt(
    id.vars = c("patient", "timepoint", "malignant_status_manual"),
    measure.vars = hla_genes,
    variable.name = "Symbol",
    value.name = "logcounts"
  ) %>%
  dplyr::mutate(timepoint = factor(timepoint),
                Symbol = factor(Symbol)) %>%
  dplyr::left_join(nonmalignant_de_hla %>%
                     dplyr::mutate(malignant_status_manual = "nonmalignant") %>%
                     dplyr::select(c(FDR, Symbol, logFC.T1))) %>%
  dplyr::left_join(malignant_de_hla %>%
                     dplyr::mutate(malignant_status_manual = "malignant") %>%
                     dplyr::select(c(FDR, Symbol, logFC.T1)))

hla_exprs_annotated_summary <- hla_exprs_annotated %>%
  dplyr::group_by(patient, timepoint, malignant_status_manual, Symbol, logFC.T1) %>%
  dplyr::summarise(logcount_mean=mean(logcounts)) %>%
  dplyr::ungroup()

hla_lineplot <- ggplot(hla_exprs_annotated, aes(x=timepoint, y=logcounts, colour=Symbol)) + 
  stat_summary(aes(colour=Symbol), fun.data = function(x) mean_se(x, mult = 10), geom = "errorbar", width = 0.05) +
  geom_line(data = hla_exprs_annotated_summary, aes(x=timepoint, y=logcount_mean,
                                                    colour=Symbol, group=Symbol)) + 
  geom_point(data = hla_exprs_annotated_summary, aes(x=timepoint, y=logcount_mean, size = -logFC.T1)) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Timepoint") + 
  ylab("Expression") + 
  scale_x_discrete(expand = c(0.2, 0.2)) + 
  guides(size = guide_legend(title = "-logFC")) + 
  scale_colour_manual(values = categorical_palettes$hla_genes) + 
  facet_wrap(~ malignant_status_manual) + 
  stripped_theme()

# Cell cycle T2 vs. T1

cycling_stats <- with(colData(sce), table(timepoint, celltype_full, Cell_Cycle)) %>% 
  as.data.frame %>%
  dplyr::rename(count=Freq) %>%
  dplyr::group_by(timepoint, celltype_full) %>% 
  dplyr::mutate(total_count=sum(count)) %>%
  dplyr::ungroup()

cycling_summary <- cycling_stats %>% 
  dplyr::filter(Cell_Cycle %in% c("S", "G2M")) %>% 
  dplyr::group_by(timepoint, celltype_full, total_count) %>%
  dplyr::summarise(cycling_count=sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cycling_prop=cycling_count/total_count) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(total_prop=total_count/sum(total_count)) %>%
  dplyr::ungroup()

min_cell_count <- 15

removed_cell_types <- (cycling_summary %>%
  dplyr::filter(total_count < min_cell_count))$celltype_full

cycling_plot <- ggplot(cycling_summary %>% dplyr::filter(!celltype_full %in% removed_cell_types),
       aes(x=timepoint, y=cycling_prop, colour=celltype_full)) +
  geom_point(aes(size=total_prop)) +
  geom_line(aes(group=celltype_full)) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() +
  xlab("Timepoint") + 
  ylab("Proportion of cells in S/G2/M") + 
  scale_size_continuous(trans = "log", breaks = c(0.01, 0.1, 0.5)) + 
  scale_colour_manual(values = categorical_palettes$celltype) + 
  guides(colour=guide_legend(title="Celltype"),
         size = guide_legend(title = "Proportion"))


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

# Build combined plot
proportion_plots <- cowplot::plot_grid(tcell_proportion_plot, bcell_proportion_plot, labels = c('b', 'c'),
                                       ncol = 1, nrow = 2)
top_row <- cowplot::plot_grid(dr_tcell, proportion_plots, labels = c('a', ''), ncol = 2, nrow = 1, rel_widths = c(0.7, 0.3))
top_row_annotated <- cowplot::plot_grid(top_row, celltype_timepoint_legend, labels = c('', ''), nrow = 2, ncol = 1,
                                        rel_heights = c(0.8, 0.2))
middle_row <- cowplot::plot_grid(cd8_activation_plot, labels = c('d'), ncol =1, nrow =1)
bottom_row <- cowplot::plot_grid(hla_lineplot, cycling_plot, rel_widths = c(0.6, 0.4), ncol = 2, nrow = 1, labels = c('e', 'f'))


final_plot <- cowplot::plot_grid(top_row_annotated, middle_row, bottom_row, 
                                 labels = c('', '', ''), 
                                 ncol = 1, 
                                 nrow = 3,
                                 rel_heights = c(0.45, 0.25, 0.3))

# Plot final plot
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")


