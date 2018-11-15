# Differential expression accounting for malignant status and timepoint (for B cells)

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

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Perform DE tests on B cells")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for B cells")
parser$add_argument('--min_gene_counts', type='integer', default = 500,
                    help="Minimum number of gene counts to filter at")
parser$add_argument('--patients', type='character', nargs ='+',
                    help="Patients to use", default = NULL)
parser$add_argument('--ngene', type='integer',
                    help="Number of top genes to use", default = 50)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for DE results.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

object <- sce %>%
  scater::mutate(malignant_status_manual = factor(malignant_status_manual, levels = c("nonmalignant", "malignant")))
object <- object[rowData(object)$total_counts >= args$min_gene_counts,]

# Filter by patient
if (!is.null(args$patients)) {
  object <- object %>%
    scater::filter(patient %in% unlist(args$patients))
}

coefficients <- c("malignant_status_manualmalignant",
                  "malignant_status_manualmalignant:timepointT2",
                  "timepointT2")

labels <- c("malignant",
            "malignant-timepoint",
            "timepoint")

de_res <- lapply(coefficients, function(coef) {
  voom_res <- de_analysis(object, de_method = "voom", formula = ~ malignant_status_manual + timepoint + malignant_status_manual*timepoint, 
                          cluster = NULL, 
                          filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = coef, 
                          min_gene_counts = args$min_gene_counts)
  
  up_genes <- (voom_res$de_table %>% dplyr::filter(logFC > 0, FDR < 0.05) %>% dplyr::arrange(-logFC))$Symbol
  down_genes <- (voom_res$de_table %>% dplyr::filter(logFC < 0, FDR < 0.05) %>% dplyr::arrange(logFC))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(object), voom_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
  
  up_genes <- up_genes %>% head(args$ngene)
  down_genes <- down_genes %>% head(args$ngene)
  
  enriched_paths <- lapply(list(up_genes, down_genes), function(genes) {
    gene_entrez <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
    paths <- enrichPathway(gene=unname(gene_entrez),pvalueCutoff=0.05, readable=TRUE, 
                           universe = unname(background_genes))
    
    return(paths)
  })
  names(enriched_paths) <- c("up", "down")
  
  return(list(gene=voom_res$de_table,
              pathway=enriched_paths))
})

names(de_res) <- labels

# Save DE result (pathway and gene tables)
saveRDS(de_res, args$outfname)

cat("Completed.\n")



