# CellBench plot

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(Matrix)
library(ggrastr)
library(mclust)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "CellBench plot")
parser$add_argument('--sce_tian_3_10x', metavar='FILE', type='character',
                    help="SCE of 10x data")
parser$add_argument('--sce_tian_3_celseq', metavar='FILE', type = 'character',
                    help="SCE of CELSeq2 data")
parser$add_argument('--sce_tian_3_dropseq', metavar='FILE', type = 'character',
                    help="SCE of Dropseq data")
parser$add_argument('--sce_tian_3_mixture', metavar='FILE', type='character',
                    help="SCE of CELSeq2 mixture")
parser$add_argument('--fit_tian_pure', metavar='FILE', type='character',
                    help="CellAssign fit to pure data")
parser$add_argument('--fit_tian_20', metavar='FILE', type='character',
                    help="CellAssign fit to mixture, with 20 marker genes")
parser$add_argument('--fit_tian_30', metavar='FILE', type='character',
                    help="CellAssign fit to mixture, with 30 marker genes")
parser$add_argument('--fit_tian_50', metavar='FILE', type='character',
                    help="CellAssign fit to mixture, with 50 marker genes")
parser$add_argument('--cell_lines', type='character', nargs = '+',
                    help="Cell lines to use")
parser$add_argument('--dimreduce_type', type='character',
                    help="Type of reduced dimension plot", choices = c("UMAP", "PCA", "TSNE"))
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

fit_tian_pure <- readRDS(args$fit_tian_pure)
fit_tian20 <- readRDS(args$fit_tian_20)
fit_tian30 <- readRDS(args$fit_tian_30)
fit_tian50 <- readRDS(args$fit_tian_50)

cell_lines <- unlist(args$cell_lines)

tian_3_10x <- readRDS(args$sce_tian_3_10x)
tian_3_celseq <- readRDS(args$sce_tian_3_celseq)
tian_3_Dropseq <- readRDS(args$sce_tian_3_dropseq)
tian_3_mixture <- readRDS(args$sce_tian_3_mixture)

colData(tian_3_10x)$sample_id <- '10x'
colData(tian_3_celseq)$sample_id <- 'CELSeq2'
colData(tian_3_Dropseq)$sample_id <- 'Dropseq'
colData(tian_3_mixture)$sample_id <- 'mixture'

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

sce_tian_pure <- merge_sces(list(tian_3_10x, tian_3_celseq, tian_3_Dropseq))
sce_tian_mix <- merge_sces(list(tian_3_celseq, tian_3_mixture))

sce_tian_pure$cellassign_class <- fit_tian_pure$cell_type

sce_tian_pure <- normalize(sce_tian_pure)
sce_tian_pure <- runPCA(sce_tian_pure, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
sce_tian_pure <- runTSNE(sce_tian_pure, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
sce_tian_pure <- runUMAP(sce_tian_pure, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)

tian_celltype <- get_colour_palette(sce_tian_pure$cell_line)

cell_line_plot <- plotReducedDim(sce_tian_pure, use_dimred = args$dimreduce_type, colour_by = "cell_line") + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  guides(fill = guide_legend(title = "Cell line", ncol = 2)) + 
  scale_fill_manual(values = tian_celltype) +
  ggtitle("Cell line") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature()

sample_plot <- plotReducedDim(sce_tian_pure, use_dimred = args$dimreduce_type, colour_by = "sample_id") + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  guides(fill = guide_legend(title = "Method", ncol = 2)) + 
  ggtitle("Method") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature()

cellassign_plot <- plotReducedDim(sce_tian_pure, use_dimred = args$dimreduce_type, colour_by = "cellassign_class") + 
  xlab(paste0(args$dimreduce_type, "-1")) + 
  ylab(paste0(args$dimreduce_type, "-2")) + 
  guides(fill = guide_legend(title = "Cellassign type", ncol = 2)) + 
  scale_fill_manual(values = tian_celltype) +
  ggtitle("CellAssign") + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature()

## SCE mix

for (cl in cell_lines) {
  colData(sce_tian_mix)[,cl][sce_tian_mix$cell_line == cl] <- 9
  colData(sce_tian_mix)[,cl][sce_tian_mix$cell_line != cl] <- 0
}

fit_objects <- list('20'=fit_tian20, '30'=fit_tian30, '50'=fit_tian50)
cellassign_probs <- plyr::rbind.fill(lapply(seq_along(fit_objects), function(i) {
  fit <- fit_objects[[i]]
  gammas <- reshape2::melt(fit$mle_params$gamma) %>%
    dplyr::rename(cell_id=Var1, cell_line=Var2, cellassign_prob=value) %>%
    dplyr::mutate(group=names(fit_objects)[i])
  return(gammas)
}))
cellassign_probs$group <- factor(cellassign_probs$group)

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
                sample_id == "mixture")

## CellBench probability plot
cellbench_prob_plot <- ggplot(prob_df, aes(x=factor(num_cells), y=cellassign_prob, fill = factor(group))) + 
  geom_boxplot(alpha = 0.7, position = "dodge", outlier.size = -1) + 
  geom_point(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0)) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  facet_wrap(~ cell_line, ncol = 1) + 
  stripped_theme(strip_face = "bold") + 
  guides(fill = guide_legend(title = "Number of marker genes")) + 
  xlab("# cells in 9-cell barcode") + 
  ylab("CellAssign probability")

entropies <- prob_df %>% 
  dplyr::group_by(cell_id, group, sample_id) %>% 
  dplyr::summarise(pure=ifelse(max(num_cells) == 9, "Pure", "Mixture"), entropy=entropy.empirical(cellassign_prob))

pvals <- compute_pvals_subsets(entropies, 
                               facet_vars = c("group"), 
                               formula = entropy ~ pure, 
                               corfun = wilcox.test, 
                               output = "p.value")

entropy_plot <- ggplot(entropies, aes(x=pure, y=entropy)) + 
  geom_boxplot(position = "dodge", outlier.size = -1) + 
  geom_point(alpha = 0.5, position = position_jitter(width = 0.2, height = 0)) + 
  theme_bw() + 
  theme_Publication() + 
  theme_nature() + 
  facet_wrap(~ group, ncol = 3) + 
  stripped_theme(strip_face = "bold") + 
  xlab("Pseudocell type") + 
  ylab("Entropy") + 
  geom_text(data = pvals, aes(label=p.adj.text), parse = TRUE, x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3) + 
  ggtitle("Number of markers")


modality_plots <- cowplot::plot_grid(cell_line_plot, 
                                     sample_plot, 
                                     cellassign_plot, 
                                     ncol = 3, 
                                     align = 'hv', 
                                     axis = 'tblr')

final_plot <- cowplot::plot_grid(modality_plots,
                                 cellbench_prob_plot,
                                 entropy_plot,
                                 nrow = 3,
                                 rel_heights = c(1, 1.2, 1))

# Plot final plot
pdf(args$outfname, width = 10, height = 12, useDingbats = FALSE)
plot(final_plot)
dev.off()

cat("Completed.\n")

