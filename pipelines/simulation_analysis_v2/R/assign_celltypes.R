
library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(SC3)
library(Matrix)
library(mclust)
library(multiClust)
library(limma)
library(edgeR)
library(infotheo)
library(future)
library(tensorflow)
library(crch)
library(splitstackshape)
library(cellassign)
library(SCINA)

library(scrna.utils)
library(cellassign)
library(cellassign.utils)
library(scrna.sceutils)
library(Seurat)

library(methods)
library(argparse)

parser <- ArgumentParser(description = "Cluster SCE with CellAssign")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Input SCE object")
parser$add_argument('--out_clusters', metavar='FILE', type = 'character',
                    help="Output directory for clustering results.")
parser$add_argument('--out_deltas', metavar='FILE', type = 'character',
                    help="Output directory for delta results.")
parser$add_argument('--out_fit', metavar='FILE', type = 'character',
                    help="Output directory for CellAssign fit")
parser$add_argument('--conda_env', type = 'character',
                    help="Conda environment", default = 'r-tensorflow')
parser$add_argument('--fc_percentile', type = 'double',
                    help="Fold change percentile cutoff", default = 0.95)
parser$add_argument('--expr_percentile', type = 'double',
                    help="Mean expression percentile cutoff", default = 0.9)
parser$add_argument('--frac_genes', type = 'double',
                    help="Fraction of genes to use", default = 1)
parser$add_argument('--test_proportion', type = 'double',
                    help="Test proportion (for supervised methods)", default = 0.5)
parser$add_argument('--max_genes', type = 'integer',
                    help="Maximum number of marker genes per cell type", default = 15)
parser$add_argument('--method', type = 'character',
                    help="Assigning method", default = "cellassign")
parser$add_argument('--wrong_prop', type = 'double',
                    help="Wrong marker proportion", default = 0)
args <- parser$parse_args()

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

fc_percentile <- args$fc_percentile
expr_percentile <- args$expr_percentile
frac_genes <- args$frac_genes
test_proportion <- args$test_proportion
max_genes <- args$max_genes
wrong_marker_proportion <- args$wrong_prop
method <- args$method

sce_total <- readRDS(args$sce)

## Split SCE into training and test sets
coldat <- as.data.frame(colData(sce_total))
sets <- splitstackshape::stratified(coldat, group = "Group", size = (1-test_proportion), bothSets = TRUE)

rowData(sce_total)$symbol <- rowData(sce_total)$Gene
rowData(sce_total)$feature_symbol <- rowData(sce_total)$Gene
rowData(sce_total)$Symbol <- rowData(sce_total)$Gene

counts(sce_total) <- as.matrix(counts(sce_total))

sce_train <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP1$Cell))]
sce <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP2$Cell))]


markers_to_use <- select_markers(sce_total, 
                                 percentile_fc = fc_percentile, 
                                 percentile_meanexpr = expr_percentile, 
                                 frac_genes = frac_genes, 
                                 max_genes_per_class = max_genes)
rho <- create_rho_matrix(sce_total, 
                         markers_to_use, 
                         wrong_proportion = wrong_marker_proportion)

s <- sizeFactors(sce)
B <- 20
X <- NULL

marker_list <- lapply(colnames(rho), function(x) {
  genes <- rownames(rho)[rho[,x] == 1]
  if (length(genes) > 0) {
    return(genes)
  }
})
names(marker_list) <- colnames(rho)
marker_list <- marker_list[sapply(marker_list, function(x) !is.null(x))]

gene_subset <- markers_to_use
n_markers <- length(markers_to_use)

if (method == "cellassign") {
  res <- cellassign(exprs_obj = sce[rownames(rho),], 
                    s = s, 
                    marker_gene_info = rho, 
                    X = X, 
                    B = B, 
                    shrinkage = TRUE, 
                    verbose = FALSE, 
                    rel_tol_em = 1e-5,
                    min_delta = 0, 
                    num_runs = 5)
  
  sce$cluster <- str_replace(res$cell_type, "^DEFac", "")
  
  delta_compare_res <- delta_compare(sce[rownames(rho),], res, colour_by = "Group", shape_by = "high_expr")
  delta_compare_df <- delta_compare_res$df
} else {
  expr_mat <- as.matrix(logcounts(sce))
  rownames(expr_mat) <- rowData(sce)$Symbol
  
  res <- SCINA(exp = expr_mat, 
               signatures = marker_list, 
               max_iter = 1000,
               rm_overlap = 0)
  
  sce$cluster <- str_replace(res$cell_labels, "^DEFac", "")
  
  delta_compare_df <- data.frame()
}

cluster_df <- data.frame(cell=as.character(colData(sce)$Cell), cluster=sce$cluster)

clusters_annotated <- data.frame(sce@metadata$param_df,
                                 clustering_method=method,
                                 frac_genes=frac_genes,
                                 max_genes=max_genes,
                                 fc_percentile=fc_percentile,
                                 expr_percentile=expr_percentile,
                                 num_markers=n_markers,
                                 test_proportion=test_proportion,
                                 cluster_df,
                                 stringsAsFactors = FALSE)

if (nrow(delta_compare_df) > 0) {
  deltas_annotated <- data.frame(sce@metadata$param_df,
                                 clustering_method=method,
                                 frac_genes=frac_genes,
                                 max_genes=max_genes,
                                 fc_percentile=fc_percentile,
                                 expr_percentile=expr_percentile,
                                 num_markers=n_markers,
                                 test_proportion=test_proportion,
                                 delta_compare_df,
                                 stringsAsFactors = FALSE)
} else {
  deltas_annotated <- delta_compare_df
}

## Write outputs

write.table(clusters_annotated, file = args$out_clusters, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(deltas_annotated, file = args$out_deltas, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
saveRDS(res, args$out_fit)


