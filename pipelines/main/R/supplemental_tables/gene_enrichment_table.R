# Supplemental table containing gene enrichment results

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(Matrix)
library(xlsx)
library(ImportExport)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create gene enrichment supplemental table")
parser$add_argument('--de_timepoint_dir', metavar='DIR', type='character',
                    help="Path to DE timepoint results")
parser$add_argument('--de_malignant_timepoint_dir', metavar = 'DIR', type = 'character',
                    help="Path to DE malignant-timepoint results")
parser$add_argument('--de_timepoint_fgsea_dir', metavar = 'DIR', type = 'character',
                    help="Path to DE fgsea dir")
parser$add_argument('--de_site_dir', metavar='DIR', type='character',
                    help="Directory to DE site results")
parser$add_argument('--de_site_fgsea_dir', metavar='DIR', type='character',
                    help="Directory to DE site fgsea results")
parser$add_argument('--de_epithelial_dir', metavar='DIR', type='character',
                    help="Directory to DE epithelial clusters results")
parser$add_argument('--padj_threshold', type = 'double',
                    help="Threshold to filter pathways at", default = 0.05)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()

## FOLLICULAR

de_timepoint_dir <- args$de_timepoint_dir
de_malignant_timepoint_dir <- args$de_malignant_timepoint_dir
de_timepoint_fgsea_dir <- args$de_timepoint_fgsea_dir
padj_threshold <- args$padj_threshold

de_timepoint_files <- Sys.glob(file.path(de_timepoint_dir, "*", "*"))
de_malignant_timepoint_files <- Sys.glob(file.path(de_malignant_timepoint_dir, "b_cells", "*"))
de_timepoint_fgsea_files <- Sys.glob(file.path(de_timepoint_fgsea_dir, "malignant", "*"))

de_timepoint_celltypes_include <- c("b", "cytotoxic", "helper", "follicular_helper", "malignant")
de_timepoint_files <- de_timepoint_files[basename(dirname(de_timepoint_files)) %in% de_timepoint_celltypes_include]

## HGSC

de_site_dir <- args$de_site_dir
de_site_fgsea_dir <- args$de_site_fgsea_dir
de_epithelial_dir <- args$de_epithelial_dir

de_site_fgsea_files <- Sys.glob(file.path(de_site_fgsea_dir, "epithelial", "*"))
de_epithelial_files <- Sys.glob(file.path(de_epithelial_dir, "*"))

read_gene_result <- function(de_files, malignant_timepoint = FALSE) {
  res <- plyr::rbind.fill(lapply(de_files, function(f) {
    de_res <- readRDS(f)
    celltype <- basename(dirname(f))
    patient <- tools::file_path_sans_ext(basename(f))
    
    if (malignant_timepoint) {
      gene_table <- de_res$malignant$gene
    } else {
      gene_table <- de_res$gene
    }
    
    enriched_genes <- gene_table %>%
      dplyr::filter(FDR <= padj_threshold) %>%
      dplyr::mutate(celltype=celltype,
                    patient=patient)
    
    return(enriched_genes)
  }))
  return(res)
}

## DE timepoint genes
timepoint_genes <- read_gene_result(de_timepoint_files)

## DE timepoint pathways for malignant cells (i.e. hallmark, fgsea)
hallmark_genes <- read_gene_result(de_timepoint_fgsea_files) %>%
  dplyr::filter(contrast == "T2") %>%
  dplyr::select(-c(contrast, logFC.T2))

## DE malignant-timepoint pathways (for only the malignant contrast term)
## NOT SHOWN AT THE MOMENT (since only part of this mentioned in the text is for genes, not pathways)
malignant_timepoint_genes <- read_gene_result(de_malignant_timepoint_files, 
                                              malignant_timepoint = TRUE)


#### HGSC section


## Hallmark for malignant cells
hgsc_hallmark_genes <- read_gene_result(de_site_fgsea_files) %>%
  dplyr::filter(contrast == "Left ovary") %>% 
  dplyr::select(-c(logFC.Left.ovary, Top))

## Hallmark for malignant cell clusters
hgsc_epithelial_cluster_hallmark_genes <- read_gene_result(de_epithelial_files) %>%
  dplyr::rename(cluster=contrast) %>%
  dplyr::select(-c(celltype))


# Write supplemental table xlsx
excel_export(list(timepoint_genes,
                  hallmark_genes,
                  hgsc_hallmark_genes,
                  hgsc_epithelial_cluster_hallmark_genes),
             file = args$outfname,
             table_names = c("T2 vs. T1 reactome genes (by celltype)",
                             "T2 vs. T1 hallmark genes (malignant B cells)",
                             "Left ovary vs. right ovary hallmark genes (epithelial cells)",
                             "Epithelial cell clusters hallmark genes"))

cat("Completed.\n")


