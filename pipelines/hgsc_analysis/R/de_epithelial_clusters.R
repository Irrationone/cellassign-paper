# Differential expression between epithelial clusters

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(gage)
library(limma)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Perform DE tests between epithelial clusters")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--cluster_col', type='character', 
                    help="Cluster column", default = "cluster")
parser$add_argument('--patients', type='character', nargs ='+',
                    help="Patients to use", default = NULL)
parser$add_argument('--gene_set_file', type='character', metavar = "FILE",
                    help="Gene set file path, if using fgsea", default = NULL)
parser$add_argument('--min_gene_counts', type='integer', default = 500,
                    help="Minimum number of gene counts to filter at")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for DE tables.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
method_gene <- args$method_gene
cluster_col <- args$cluster_col

sce_filtered <- sce[,!is.na(colData(sce)[,cluster_col])]

# Filter by patient
if (!is.null(args$patients)) {
  sce_filtered <- sce_filtered %>%
    scater::filter(patient %in% unlist(args$patients))
}

scran_res <- de_analysis(sce_filtered, de_method = "scran", formula = NULL, cluster = colData(sce_filtered)[,cluster_col], 
                         filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = NULL, 
                         min_gene_counts = args$min_gene_counts, full.stats = TRUE)
de_table <- scran_res$de_table

if (!is.null(args$gene_set_file)) {
  pathway_list <- gmtPathways(args$gene_set_file)
} else {
  stop("Must specify a gene set file.")
}

clust1 <- as.character(unique(de_table$contrast))
pathway_combined <- plyr::rbind.fill(lapply(clust1, function(x) {
  de_subset <- de_table %>% 
    dplyr::filter(contrast == x)
  
  other_clusts <- setdiff(clust1, x)
  pathway_df <- plyr::rbind.fill(lapply(other_clusts, function(y) {
    logfc_col <- paste0("logFC", ".", y)
    df <- data.frame(Symbol=de_subset$Symbol, logFC=de_subset[,logfc_col])
    df <- df %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarize(logFC=mean(logFC)) %>%
      dplyr::ungroup()
    
    fgsea_res <- fgsea(pathways=pathway_list, 
                       stats = deframe(df), 
                       nperm = 10000)
    
    return(data.frame(clust1=x, clust2=y, fgsea_res))
  }))
  return(pathway_df)
}))


# Save DE result (pathway and gene tables)
saveRDS(list(gene=de_table, pathway=pathway_combined), args$outfname)

cat("Completed.\n")



