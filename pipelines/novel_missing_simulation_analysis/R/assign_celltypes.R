
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
parser$add_argument('--max_genes', type = 'integer',
                    help="Maximum number of marker genes per cell type", default = 15)
parser$add_argument('--method', type = 'character',
                    help="Assigning method", default = "cellassign")
parser$add_argument('--wrong_prop', type = 'double',
                    help="Wrong marker proportion", default = 0)
parser$add_argument('--data_subset', type = 'character', nargs ='+',
                    help="Celltypes to use in the data")
parser$add_argument('--marker_subset', type = 'character', nargs ='+',
                    help="Celltypes to search for")
args <- parser$parse_args()

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

fc_percentile <- args$fc_percentile
expr_percentile <- args$expr_percentile
frac_genes <- args$frac_genes
test_proportion <- args$test_proportion
max_genes <- args$max_genes
wrong_marker_proportion <- args$wrong_prop
method <- args$method
data_subset_types <- unlist(args$data_subset)
marker_subset_types <- unlist(args$marker_subset)

sce <- readRDS(args$sce)

rowData(sce)$symbol <- rowData(sce)$Gene
rowData(sce)$feature_symbol <- rowData(sce)$Gene
rowData(sce)$Symbol <- rowData(sce)$Gene

all_groups <- sort(unique(sce$Group))

counts(sce) <- as.matrix(counts(sce))

markers_to_use <- select_markers(sce, 
                                 percentile_fc = fc_percentile, 
                                 percentile_meanexpr = expr_percentile, 
                                 frac_genes = frac_genes, 
                                 max_genes_per_class = max_genes)
rho <- create_rho_matrix(sce, 
                         markers_to_use, 
                         wrong_proportion = wrong_marker_proportion)
colnames(rho) <- str_replace(colnames(rho), "^DEFac", "")

# Subset data
sce <- sce %>%
  scater::filter(Group %in% data_subset_types)

# Subset rho
rho <- rho[,marker_subset_types]

marker_list <- lapply(colnames(rho), function(x) {
  genes <- rownames(rho)[rho[,x] == 1]
  if (length(genes) > 0) {
    return(genes)
  }
})
names(marker_list) <- colnames(rho)

print(marker_list)

gene_subset <- markers_to_use
n_markers <- length(markers_to_use)

s <- sizeFactors(sce)
B <- 20
X <- NULL

if (method == "cellassign") {
  # Add other group to rho
  rho_exp <- marker_list_to_mat(marker_list, include_other = TRUE)
  
  res <- cellassign(exprs_obj = sce[rownames(rho_exp),], 
                    s = s, 
                    marker_gene_info = rho_exp, 
                    X = X, 
                    B = B, 
                    shrinkage = TRUE, 
                    verbose = FALSE, 
                    rel_tol_em = 1e-5,
                    min_delta = 0, 
                    num_runs = 5)
  
  sce$cluster <- plyr::mapvalues(res$cell_type, from = c("other"), to = c("unknown"))
} else {
  expr_mat <- as.matrix(logcounts(sce))
  rownames(expr_mat) <- rowData(sce)$Symbol
  
  res <- SCINA(exp = expr_mat, 
               signatures = marker_list, 
               max_iter = 1000,
               rm_overlap = 0, 
               allow_unknown = 1)
  
  sce$cluster <- res$cell_labels
}

cluster_df <- data.frame(cell=as.character(colData(sce)$Cell), cluster=sce$cluster)

clusters_annotated <- data.frame(sce@metadata$param_df,
                                 clustering_method=method,
                                 frac_genes=frac_genes,
                                 max_genes=max_genes,
                                 fc_percentile=fc_percentile,
                                 expr_percentile=expr_percentile,
                                 num_markers=n_markers,
                                 wrong_marker_proportion=wrong_marker_proportion,
                                 all_types=paste(all_groups, collapse = ","),
                                 data_subset_types=paste(data_subset_types, collapse = ","),
                                 marker_subset_types=paste(marker_subset_types, collapse = ","),
                                 n_data_types=length(data_subset_types),
                                 n_marker_types=length(marker_subset_types),
                                 data_subset_markers=all(data_subset_types %in% marker_subset_types),
                                 marker_subset_data=all(marker_subset_types %in% data_subset_types),
                                 cluster_df,
                                 stringsAsFactors = FALSE)

## Write outputs

write.table(clusters_annotated, file = args$out_clusters, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
saveRDS(res, args$out_fit)


