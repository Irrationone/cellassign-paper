# Supplemental table containing pathway enrichment results

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

parser <- ArgumentParser(description = "Create pathway enrichment supplemental table")
parser$add_argument('--de_timepoint_dir', metavar='DIR', type='character',
                    help="Path to DE timepoint results")
parser$add_argument('--de_malignant_timepoint_dir', metavar = 'DIR', type = 'character',
                    help="Path to DE malignant-timepoint results")
parser$add_argument('--de_timepoint_fgsea_dir', metavar = 'DIR', type = 'character',
                    help="Path to DE fgsea dir")
parser$add_argument('--de_site', metavar='DIR', type='character',
                    help="Directory to DE site results")
parser$add_argument('--de_site_fgsea', metavar='DIR', type='character',
                    help="Directory to DE site fgsea results")
parser$add_argument('--de_epithelial', metavar='DIR', type='character',
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

de_site_dir <- args$de_site
de_site_fgsea_dir <- args$de_site_fgsea
de_epithelial_dir <- args$de_epithelial

de_site_fgsea_files <- Sys.glob(file.path(de_site_fgsea_dir, "epithelial", "*"))
de_epithelial_files <- Sys.glob(file.path(de_epithelial_dir, "*"))

read_pathway_result <- function(de_files) {
  res <- plyr::rbind.fill(lapply(de_files, function(f) {
    de_res <- readRDS(f)
    celltype <- basename(dirname(f))
    patient <- tools::file_path_sans_ext(basename(f))
    
    enriched_pathways <- plyr::rbind.fill(lapply(seq_along(de_res$pathway), function(i) {
      direction <- names(de_res$pathway)[i]
      df <- (de_res$pathway[[i]] %>%
               filter_pathway_result())@result %>%
        dplyr::filter(p.adjust <= padj_threshold)
      
      df <- df %>% 
        dplyr::mutate(direction = direction)
      
      return(df)
    }))
    
    enriched_pathways <- enriched_pathways %>%
      dplyr::mutate(celltype=celltype,
                    patient=patient)
    
    return(enriched_pathways)
  }))
  return(res)
}

## DE timepoint pathways
timepoint_pathways <- read_pathway_result(de_timepoint_files)

## DE timepoint pathways for malignant cells (i.e. hallmark, fgsea)
hallmark_pathways <- plyr::rbind.fill(lapply(de_timepoint_fgsea_files, function(f) {
  de_res <- readRDS(f)
  celltype <- basename(dirname(f))
  patient <- tools::file_path_sans_ext(basename(f))
  
  enriched_pathways <- de_res$pathway %>%
    dplyr::mutate(celltype=celltype,
                  patient=patient) %>%
    dplyr::filter(padj <= padj_threshold,
                  size >= 2)
  
  return(enriched_pathways)
}))
hallmark_pathways$leadingEdge <- sapply(hallmark_pathways$leadingEdge,
                                        function(x) paste(x, collapse = ", "))

## DE malignant-timepoint pathways (for only the malignant contrast term)
## NOT SHOWN AT THE MOMENT (since only part of this mentioned in the text is for genes, not pathways)
malignant_timepoint_pathways <- plyr::rbind.fill(lapply(de_malignant_timepoint_files, function(f) {
  de_res <- readRDS(f)
  patient <- tools::file_path_sans_ext(basename(f))
  
  enriched_pathways <- plyr::rbind.fill(lapply(seq_along(de_res$malignant$pathway), function(i) {
    direction <- names(de_res$malignant$pathway)[i]
    df <- (de_res$malignant$pathway[[i]] %>%
             filter_pathway_result())@result %>%
      dplyr::filter(p.adjust <= padj_threshold)
    
    df <- df %>% 
      dplyr::mutate(direction = direction)
    
    return(df)
  }))
  
  enriched_pathways <- enriched_pathways %>%
    dplyr::mutate(patient=patient)
  
  return(enriched_pathways)
}))



#### HGSC section


## Hallmark for malignant cells
hgsc_hallmark_pathways <- plyr::rbind.fill(lapply(de_site_fgsea_files, function(f) {
  de_res <- readRDS(f)
  celltype <- basename(dirname(f))
  patient <- tools::file_path_sans_ext(basename(f))
  
  enriched_pathways <- de_res$pathway %>%
    dplyr::mutate(celltype=celltype,
                  patient=patient,
                  ES = -1*ES,
                  NES = -1*NES) %>%
    dplyr::filter(padj <= padj_threshold,
                  size >= 2)
  
  return(enriched_pathways)
}))
hgsc_hallmark_pathways$leadingEdge <- sapply(hgsc_hallmark_pathways$leadingEdge,
                                             function(x) paste(x, collapse = ", "))


## Hallmark for malignant cell clusters
hgsc_epithelial_cluster_hallmark_pathways <- plyr::rbind.fill(lapply(de_epithelial_files, function(f) {
  de_res <- readRDS(f)
  patient <- tools::file_path_sans_ext(basename(f))
  
  enriched_pathways <- de_res$pathway %>%
    dplyr::mutate(patient=patient) %>%
    dplyr::filter(padj <= padj_threshold,
                  size >= 2) %>%
    dplyr::rename(cluster_to=clust1,
                  cluster_from=clust2)
  
  return(enriched_pathways)
}))
hgsc_epithelial_cluster_hallmark_pathways$leadingEdge <- sapply(hgsc_epithelial_cluster_hallmark_pathways$leadingEdge,
                                                                function(x) paste(x, collapse = ", "))



# Write supplemental table xlsx
excel_export(list(timepoint_pathways,
                  hallmark_pathways,
                  hgsc_hallmark_pathways,
                  hgsc_epithelial_cluster_hallmark_pathways),
             file = args$outfname,
             table_names = c("T2 vs. T1 reactome pathways (by celltype)",
                             "T2 vs. T1 hallmark pathways (malignant B cells)",
                             "Left ovary vs. right ovary hallmark pathways (epithelial cells)",
                             "Epithelial cell clusters hallmark pathways"))

cat("Completed.\n")


