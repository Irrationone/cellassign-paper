# Overview figure showing FL cell types

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create overview figure for FL")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--sce_raw', metavar='FILE', type='character',
                    help="Path to raw SingleCellExperiment RDS")
parser$add_argument('--patient_progression', metavar='FILE', type='character',
                    help="Patient events")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of dimensionality reduction to plot", default = "UMAP")
parser$add_argument('--winsorized_expression_threshold', type='double',
                    help="Winsorized expression threshold", default = NULL)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

sce_path <- args$sce
sce_raw_path <- args$sce_raw
sce <- readRDS(sce_path)
sce_raw <- readRDS(sce_raw_path)

cell_counts <- data.frame(raw_cellcount=sce_raw$dataset %>% table)
colnames(cell_counts) <- c("event", "cellcount")

categorical_palettes <- cat_palettes()

patient_progression <- fread(args$patient_progression, sep = "\t")
patient_progression <- patient_progression %>%
  dplyr::mutate(patient = factor(patient),
                label = ifelse(timepoint %in% c("L", "D"), "", plot_label)) %>%
  dplyr::mutate(patient = factor(patient, levels = rev(levels(patient)))) %>%
  dplyr::left_join(cell_counts)

dead_patients <- (patient_progression %>% 
                    dplyr::filter(timepoint == "D"))$patient
live_patients <- (patient_progression %>% 
                    dplyr::filter(timepoint == "L"))$patient
death_boxes <- patient_progression %>%
  dplyr::filter(timepoint == "D") %>%
  dplyr::mutate(ycenter=as.numeric(patient),
                ymin=ycenter-0.15,
                ymax=ycenter+0.15)

timepoint_plot <- ggplot(patient_progression, aes(x=years, y=patient)) + 
  geom_text(data = patient_progression %>% dplyr::filter(observed == 1),
                  aes(label=label), size = 2.5, nudge_y = 0.2) + 
  geom_text(data = patient_progression %>% dplyr::filter(observed == 1),
            aes(label=cellcount), size = 2.5, nudge_y = -0.2) +
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  xlab("Time since first biopsy (years)") + 
  ylab("") + 
  scale_fill_manual(values = categorical_palettes$dataset) +
  scale_colour_manual(values = categorical_palettes$dataset, limits = c('FL1018T0', 
                                                                        'FL1018T1', 
                                                                        'FL1018T2', 
                                                                        'FL2001T1',
                                                                        'FL2001T2',
                                                                        'FL2001D')) + 
  guides(colour = FALSE) + 
  guides(fill = FALSE) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = 'bold')) + 
  geom_line(data = patient_progression %>% dplyr::filter(patient %in% dead_patients),
            arrow = arrow(length=unit(0,"cm"), ends="last", type = "closed")) +
  geom_line(data = patient_progression %>% dplyr::filter(patient %in% live_patients),
            arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
  geom_rect(data = death_boxes,
            aes(xmin=years-0.05,xmax=years+0.05,ymin=ymin,ymax=ymax, fill=event)) + 
  geom_point(data = patient_progression %>% dplyr::filter(observed == 1 & timepoint != "D"), 
             aes(colour=event), size = 5) 

# Plot of timepoint
dr_timepoint <- plotReducedDim(sce, use_dimred = "UMAP", colour_by = "dataset", point_alpha = 0.1, add_ticks = FALSE)
dr_timepoint <- dr_timepoint + 
  guides(colour = guide_legend(title = "Sample", override.aes = list(alpha = 1,
                                                                     size = 2),
                               title.position = "top", title.hjust = 0.5),
         fill = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$dataset)
dr_timepoint$layers[[1]]$aes_params$colour <- NULL
dr_timepoint$layers[[1]]$aes_params$shape <- 16
dr_timepoint$layers[[1]]$mapping$colour <- dr_timepoint$layers[[1]]$mapping$fill

# Plot of celltype assignments
nonother_types <- sort(setdiff(unique(sce$celltype), "other"))
dr_celltype <- plotReducedDim(sce %>%
                                scater::mutate(celltype=factor(plyr::mapvalues(celltype, from = c("other"),
                                                                                    to = c("Unassigned")),
                                                                    levels = c(nonother_types, "Unassigned"))), 
                              use_dimred = "UMAP",
                              colour_by = "celltype",
                              point_alpha = 0.1, 
                              add_ticks = FALSE)
dr_celltype <- dr_celltype + 
  guides(colour = guide_legend(title = "Predicted celltype", override.aes = list(alpha = 1,
                                                                                 size = 2),
                               title.position = "top", title.hjust = 0.5),
         fill = FALSE) + 
  xlab("UMAP-1") + 
  ylab("UMAP-2") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  scale_colour_manual(values = categorical_palettes$celltype)
dr_celltype$layers[[1]]$aes_params$colour <- NULL
dr_celltype$layers[[1]]$aes_params$shape <- 16
dr_celltype$layers[[1]]$mapping$colour <- dr_celltype$layers[[1]]$mapping$fill

# Plots of marker gene expression
marker_genes <- c("CD79A", "CD3D", "CCL5", "IL7R")
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
                      use_dimred = "UMAP",
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
    xlab("UMAP-1") + 
    ylab("UMAP-2") + 
    theme_bw() + 
    theme_Publication() + 
    theme_nature() +
    scale_colour_gradientn(colours = gradient_colours, 
                           limits = expr_limits) + 
    ggtitle(mgene)
  return(p)
})

# Legends

## Expression values

marker_legend <- cellassign.utils::ggsimplelegend(expr_limits,
                                                  colour_mapping = gradient_colours,
                                                  legend_title = "Expression",
                                                  type = "continuous") + 
  theme(legend.key.width = unit(2, "lines"))
marker_legend <- cellassign.utils::extract_legend(marker_legend)

# Assemble plot

## DR plots
dr_plots <- cowplot::plot_grid(dr_timepoint, dr_celltype, ncol = 2, nrow = 1, labels = c('b', 'c'))
marker_gene_plots <- cowplot::plot_grid(plotlist = marker_plots, ncol = 4, nrow = 1, labels = c('d', '', '', ''))

final_plot <- cowplot::plot_grid(timepoint_plot, 
                                 dr_plots, 
                                 marker_gene_plots, 
                                 marker_legend,
                                 labels = c('a', '', '', ''), 
                                 ncol = 1, 
                                 nrow = 4,
                                 rel_heights = c(0.5, 1, 0.5, 0.1))

# Plot to output file
pdf(args$outfname, width = 10, height = 10, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")



