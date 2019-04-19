
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
library(SCINA)

library(scrna.utils)
library(cellassign)
library(cellassign.utils)
library(scrna.sceutils)
library(Seurat)

library(methods)
library(argparse)

parser <- ArgumentParser(description = "Cluster SCE.")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Input SCE object")
parser$add_argument('--outdir', metavar='DIR', type = 'character',
                    help="Output directory for clustering results.")
parser$add_argument('--dimreduce_method', type = 'character',
                    help="Dimensionality reduction method.")
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
parser$add_argument('--clust_method', type = 'character',
                    help="Clustering method")
parser$add_argument('--max_genes', type = 'integer',
                    help="Maximum number of marker genes per cell type", default = 15)
parser$add_argument('--wrong_prop', type = 'double',
                    help="Wrong marker proportion", default = 0)
args <- parser$parse_args()

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

## Helper functions
cluster_wrapper <- function(sce, gene_subset = NULL, dimreduce_method, clustering_method, conda_env = "r-tensorflow", conda_binary = "/home/rstudio/miniconda/bin/conda", seurat_resolution = 0.6,
                            object2 = NULL, object2_cluster_label = "Group") {
  sce_original <- sce
  if (!is.null(gene_subset)) {
    sce <- sce[gene_subset,]
    sce <- runPCA(sce, ntop = 500, ncomponents = 50, exprs_values = "logcounts")
    
    if (!is.null(object2)) {
      object2 <- object2[gene_subset,]
    }
  }
  
  if (clustering_method == "phenograph") {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, conda_env = conda_env, 
                         conda_binary = conda_binary,
                         phenograph_module = "scrnatools.methods.clustering.phenograph_analysis")
  } else if (clustering_method == "dbscan") {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, dbscan_epsilon = 0.5, dbscan_minPoints = 5)
  } else if (str_detect(clustering_method, "^seurat")) {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method,
                         seurat_resolution = seurat_resolution)
  } else if (clustering_method %in% c("Zheng_cor", "scmap")) {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method,
                         object2 = object2, object2_cluster_label = object2_cluster_label)
  } else {
    sce <- cluster_cells(sce, method = clustering_method, dimreduce_type = dimreduce_method, ncores = 1)
  }
  
  sce_original$cluster <- sce$cluster
  return(sce_original)
}

map_and_evaluate_clusters <- function(sce_sim, min_count = 0) {
  # Remove genes that are lowly expressed, and cells with no UMIs
  sce_sim <- sce_sim[rowSums(counts(sce_sim)) >= min_count,colSums(counts(sce_sim)) > 0]
  
  nclusts <- length(unique(sce_sim$cluster))
  
  if (nclusts > 1 & !all(is.na(sce_sim$cluster))) {
    sce_de <- map_clusters(sce_sim, method = "de2", FDR_cutoff = 0.05)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0)
  } else {
    max_group <- names(which.max(table(sce_sim$Group)))[1]
    
    sce_de <- sce %>% scater::mutate(inferred_group=max_group)
    sce_correlation <- map_clusters(sce_sim, method = "correlation", min_correlation = 0)
  }
  
  de_evaluation_measures <- compute_evaluation_measures(sce_de, 
                                                        truth_labels = sce_de$Group,
                                                        inferred_labels = sce_de$inferred_group)
  corr_evaluation_measures <- compute_evaluation_measures(sce_correlation, 
                                                          truth_labels = sce_correlation$Group,
                                                          inferred_labels = sce_correlation$inferred_group)
  
  evaluation_measures <- rbind(data.frame(mapping_type = "de", de_evaluation_measures, stringsAsFactors = FALSE), data.frame(mapping_type = "correlation", corr_evaluation_measures, stringsAsFactors = FALSE))
  
  return(evaluation_measures)
}


dimreduce_method <- args$dimreduce_method
fc_percentile <- args$fc_percentile
expr_percentile <- args$expr_percentile
frac_genes <- args$frac_genes
test_proportion <- args$test_proportion
clustering_method <- args$clust_method
max_genes <- args$max_genes
wrong_marker_proportion <- args$wrong_prop

## Split SCE into training and test sets
coldat <- as.data.frame(colData(sce_total))
sets <- splitstackshape::stratified(coldat, group = "Group", size = (1-test_proportion), bothSets = TRUE)

rowData(sce_total)$symbol <- rowData(sce_total)$Gene
rowData(sce_total)$feature_symbol <- rowData(sce_total)$Gene
rowData(sce_total)$Symbol <- rowData(sce_total)$Gene

counts(sce_total) <- as.matrix(counts(sce_total))

sce_train <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP1$Cell))]
sce <- sce_total[,which(as.character(colData(sce_total)$Cell) %in% as.character(sets$SAMP2$Cell))]


## Clustering code
markers_to_use <- select_markers(sce_total, 
                                 percentile_fc = fc_percentile, 
                                 percentile_meanexpr = expr_percentile, 
                                 frac_genes = frac_genes, 
                                 max_genes_per_class = max_genes)

rho <- create_rho_matrix(sce_total, markers_to_use, wrong_proportion = wrong_marker_proportion)

s <- sizeFactors(sce)
B <- 20
X <- NULL

if (clustering_method == "cellassign") {
  res <- cellassign(exprs_obj = sce[rownames(rho),], 
                    s = s, 
                    rho = rho, 
                    X = X, 
                    B = B, 
                    shrinkage = TRUE, 
                    verbose = FALSE, 
                    rel_tol_em = 1e-5, 
                    min_delta = 0, num_runs = 5)
  
  sce$cluster <- str_replace(res$cell_type, "^DEFac", "")
  
  evaluation_measures <- compute_evaluation_measures(sce, 
                                                     truth_labels = sce$Group,
                                                     inferred_labels = sce$cluster)
  
  delta_compare_res <- delta_compare(sce[rownames(rho),], res, colour_by = "Group", shape_by = "high_expr")
  delta_compare_df <- delta_compare_res$df
  cluster_df <- data.frame(cluster=sce$cluster)
} else if (clustering_method == "scina") {
  ## TODO: Implement SCINA-based clustering
  cluster_df <- data.frame(cluster=sce$cluster)
} else if (str_detect(clustering_method, "^seurat_")) {
  seurat_resolution <- as.numeric(str_replace(clustering_method, "^seurat_", ""))
  
  sce_sim_full <- cluster_wrapper(sce, gene_subset = NULL, dimreduce_method = dimreduce_method, clustering_method = "seurat", seurat_resolution = seurat_resolution)
  
  print(with(colData(sce_sim_full), table(cluster, Group)))
  evaluation_measures_full <- map_and_evaluate_clusters(sce_sim_full)
  
  sce_sim_markers <- cluster_wrapper(sce, gene_subset = markers_to_use, dimreduce_method = dimreduce_method, clustering_method = "seurat", seurat_resolution = seurat_resolution)
  
  print(with(colData(sce_sim_markers), table(cluster, Group)))
  evaluation_measures_markers <- map_and_evaluate_clusters(sce_sim_markers)
  
  evaluation_measures <- rbind(data.frame(evaluation_measures_full, gene_set="full"),
                               data.frame(evaluation_measures_markers, gene_set="markers"))
  
  sce <- sce_sim_full
  delta_compare_df <- data.frame()
  cluster_df <- data.frame(full_cluster=sce_sim_full$cluster,
                           marker_cluster=sce_sim_markers$cluster)
} else if (clustering_method %in% c("Zheng_cor", "scmap")) {
  sce_sim_full <- cluster_wrapper(sce, gene_subset = NULL, dimreduce_method = dimreduce_method, clustering_method = clustering_method, object2 = sce_train, object2_cluster_label = "Group")
  
  print(with(colData(sce_sim_full), table(cluster, Group)))
  evaluation_measures_full <- map_and_evaluate_clusters(sce_sim_full)
  
  sce_sim_markers <- cluster_wrapper(sce, gene_subset = markers_to_use, dimreduce_method = dimreduce_method, clustering_method = clustering_method, object2 = sce_train, object2_cluster_label = "Group")
  
  print(with(colData(sce_sim_markers), table(cluster, Group)))
  evaluation_measures_markers <- map_and_evaluate_clusters(sce_sim_markers)
  
  evaluation_measures <- rbind(data.frame(evaluation_measures_full, gene_set="full"),
                               data.frame(evaluation_measures_markers, gene_set="markers"))
  
  sce <- sce_sim_full
  delta_compare_df <- data.frame()
  cluster_df <- data.frame(full_cluster=sce_sim_full$cluster,
                           marker_cluster=sce_sim_markers$cluster)
} else {
  sce_sim_full <- cluster_wrapper(sce, gene_subset = NULL, dimreduce_method = dimreduce_method, clustering_method = clustering_method)
  
  print(with(colData(sce_sim_full), table(cluster, Group)))
  evaluation_measures_full <- map_and_evaluate_clusters(sce_sim_full)
  
  gc()
  sce_sim_markers <- cluster_wrapper(sce, gene_subset = markers_to_use, dimreduce_method = dimreduce_method, clustering_method = clustering_method)
  
  print(with(colData(sce_sim_markers), table(cluster, Group)))
  evaluation_measures_markers <- map_and_evaluate_clusters(sce_sim_markers)
  
  evaluation_measures <- rbind(data.frame(evaluation_measures_full, gene_set="full"),
                               data.frame(evaluation_measures_markers, gene_set="markers"))
  
  sce <- sce_sim_full
  delta_compare_df <- data.frame()
  cluster_df <- data.frame(full_cluster=sce_sim_full$cluster,
                           marker_cluster=sce_sim_markers$cluster)
}

sce@metadata$param_df <- data.frame(sce@metadata$param_df,
                                    dimreduce_method=dimreduce_method,
                                    clustering_method=clustering_method,
                                    frac_genes=frac_genes,
                                    max_genes=max_genes,
                                    fc_percentile=fc_percentile,
                                    expr_percentile=expr_percentile,
                                    num_markers=length(markers_to_use),
                                    test_proportion=test_proportion,
                                    wrong_marker_proportion=wrong_marker_proportion,
                                    stringsAsFactors = FALSE)

## Write outputs

if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

eval_measures_file <- file.path(args$outdir, "eval_measures.tsv")
param_file <- file.path(args$outdir, "params.tsv")
delta_compare_file <- file.path(args$outdir, "delta_compare.tsv")
cluster_file <- file.path(args$outdir, "clusters.tsv")

write.table(evaluation_measures, file = eval_measures_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(sce@metadata$param_df, file = param_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(delta_compare_df, file = delta_compare_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cluster_df, file = cluster_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



