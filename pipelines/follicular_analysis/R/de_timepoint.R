# Differential expression (grouped by timepoint)

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

parser <- ArgumentParser(description = "Perform DE tests on an SCE for a particular celltype(s)")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--celltypes', type='character', nargs='+',
                    help="Celltypes to include", default = NULL)
parser$add_argument('--patients', type='character', nargs ='+',
                    help="Patients to use", default = NULL)
parser$add_argument('--method_gene', type='character',
                    help="DE method to use", default = "voom")
parser$add_argument('--method_pathway', type='character',
                    help="DE method to use", default = "ReactomePA")
parser$add_argument('--gene_set_file', type='character', metavar = "FILE",
                    help="Gene set file path, if using fgsea", default = NULL)
parser$add_argument('--ngene', type='integer',
                    help="Number of top genes to use", default = 50)
parser$add_argument('--min_gene_counts', type='integer', default = 500,
                    help="Minimum number of gene counts to filter at")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for DE tables.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
method_gene <- args$method_gene

sce_filtered <- sce %>%
  scater::mutate(malignant_status_manual = factor(malignant_status_manual, levels = c("nonmalignant", "malignant")))

# Filter by celltype
if (!is.null(args$celltypes)) {
  sce_filtered <- sce_filtered %>% 
    scater::filter(celltype_full %in% unlist(args$celltypes))
}

# Filter by patient
if (!is.null(args$patients)) {
  sce_filtered <- sce_filtered %>%
    scater::filter(patient %in% unlist(args$patients))
}

cont <- "T2"
coef <- "timepointT2"
formula <- ~timepoint
cluster_col <- "timepoint"


if (method_gene == "voom") {
  voom_res <- de_analysis(sce_filtered, de_method = "voom", formula, cluster = NULL, 
                          filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = coef, 
                          min_gene_counts = args$min_gene_counts)
  de_table <- voom_res$de_table
  
  up_genes <- (voom_res$de_table %>% dplyr::filter(logFC > 0, FDR < 0.05) %>% dplyr::arrange(-logFC))$Symbol
  down_genes <- (voom_res$de_table %>% dplyr::filter(logFC < 0, FDR < 0.05) %>% dplyr::arrange(logFC))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(sce), voom_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
} else if (method_gene == "scran") {
  scran_res <- de_analysis(sce_filtered, de_method = "scran", formula = NULL, cluster = colData(sce_filtered)[,cluster_col], 
                           filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = NULL, 
                           min_gene_counts = args$min_gene_counts)
  de_table <- scran_res$de_table
  
  up_genes <- (scran_res$de_table %>% dplyr::filter(contrast == cont, FDR < 0.05, logFC.T1 > 0) %>% dplyr::arrange(-logFC.T1))$Symbol
  down_genes <- (scran_res$de_table %>% dplyr::filter(contrast == cont, FDR < 0.05, logFC.T1 < 0) %>% dplyr::arrange(logFC.T1))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(sce), scran_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
} else {
  stop("Unrecognized DE method.")
}



if (args$method_pathway == "ReactomePA") {
  up_genes <- up_genes %>% head(args$ngene)
  down_genes <- down_genes %>% head(args$ngene)
  
  enriched_paths <- lapply(list(up_genes, down_genes), function(genes) {
    gene_entrez <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
    paths <- enrichPathway(gene=unname(gene_entrez),pvalueCutoff=0.05, readable=TRUE, 
                           universe = unname(background_genes))
    
    return(paths)
  })
  names(enriched_paths) <- c("up", "down")
} else if (args$method_pathway == "fgsea") {
  # Only allow use of scran with this right now -- can use scran but need to rename columns accordingly
  #stopifnot(method_gene == "voom")
  stopifnot(method_gene == "scran")
  
  if (!is.null(args$gene_set_file)) {
    pathway_list <- gmtPathways(args$gene_set_file)
  } else {
    stop("Must specify a gene set file.")
  }
  
  de_table_summarized <- de_table %>%
    dplyr::filter(contrast == "T2") %>%
    dplyr::select(Symbol, logFC.T1) %>%
    distinct() %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(logFC = mean(logFC.T1))
  
  # de_table_summarized <- de_table %>% 
  #   dplyr::select(Symbol, t) %>% 
  #   na.omit() %>% 
  #   distinct() %>% 
  #   dplyr::group_by(Symbol) %>% 
  #   dplyr::summarise(t = mean(t))
  # 
  fgsea_res <- fgsea(pathways=pathway_list, 
                     stats = deframe(de_table_summarized), 
                     nperm = 10000)
  
  enriched_paths <- fgsea_res %>%
    as_tibble() %>%
    dplyr::arrange(desc(NES))
} else {
  stop("Unrecognized pathway enrichment method.")
}


# Save DE result (pathway and gene tables)
saveRDS(list(gene=de_table, pathway=enriched_paths), args$outfname)

cat("Completed.\n")



