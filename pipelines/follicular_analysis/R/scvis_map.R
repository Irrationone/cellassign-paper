# Map scvis

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(scvis)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Determine malignant status of B cells in FL data")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to follicular SingleCellExperiment RDS (scvis trained)")
parser$add_argument('--sce_rln', metavar='FILE', type='character',
                    help="Path to RLN SingleCellExperiment RDS")
parser$add_argument('--conda_env', type='character',
                    help="Conda environment with scanorama", default = "r-tensorflow")
parser$add_argument('--scvis_config_path', metavar='FILE', type='character',
                    help="Path to scvis config")
parser$add_argument('--scvis_output_dir', metavar='DIR', type='character',
                    help="Path to scvis temporary output directory")
parser$add_argument('--out_rln', type = 'character', metavar = 'FILE',
                    help="Output path for mapped RLN SCE.")
parser$add_argument('--out_merged', type = 'character', metavar = 'FILE',
                    help="Output path for mapped merged SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
sce_rln <- readRDS(args$sce_rln)
scvis_config_path <- args$scvis_config_path

Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

pca_res <- sce@metadata$pca_res

pca_RLN <- pca2_project(sce_rln, pca = pca_res, exprs_values = "logcounts")
reducedDim(sce_rln, "PCA2_project") <- pca_RLN

scvis_map_RLN <- scvis_map(sce_rln, 
                           output_dir = args$scvis_output_dir, 
                           use_reducedDim = TRUE,
                           reducedDim_name = "PCA2_project", 
                           scvis_config_path)

merge_sces_scvis <- function(sce1, sce2, dimreduce1, dimreduce2, name1 = "data1", name2 = "data2", common_name = "scvis", n_dims = 2, merge_rowdat = TRUE, merge_reducedDim = TRUE, keep_cols = c()) {
  if (merge_rowdat) {
    rowData(sce1) <- subset(rowData(sce1), select = c("ID", "Symbol"))
    rowData(sce2) <- subset(rowData(sce2), select = c("ID", "Symbol"))
  }
  
  for (col in keep_cols) {
    if (!col %in% colnames(colData(sce1))) {
      colData(sce1)[,col] <- NA
    }
    
    if (!col %in% colnames(colData(sce2))) {
      colData(sce2)[,col] <- NA
    }
  }
  
  common_cols <- intersect(colnames(colData(sce1)), colnames(colData(sce2)))
  colData(sce1) <- subset(colData(sce1), select = common_cols)
  colData(sce2) <- subset(colData(sce2), select = common_cols)
  colData(sce1)[,"set_name"] <- name1
  colData(sce2)[,"set_name"] <- name2
  
  common_assays <- intersect(names(assays(sce1)), names(assays(sce2)))
  assays(sce1) <- assays(sce1)[common_assays]
  assays(sce2) <- assays(sce2)[common_assays]
  
  if (merge_reducedDim) {
    reducedDims(sce1) <- SimpleList(reducedDim(sce1, dimreduce1)[,1:n_dims])
    reducedDims(sce2) <- SimpleList(reducedDim(sce2, dimreduce2)[,1:n_dims])
    names(reducedDims(sce1)) <- common_name
    names(reducedDims(sce2)) <- common_name
  } else {
    reducedDims(sce1) <- NULL
    reducedDims(sce2) <- NULL
  }
  
  sce_merged <- cbind(sce1, sce2)
  return(sce_merged)
}

scvis_combined <- merge_sces_scvis(sce, 
                                   scvis_map_RLN, 
                                   dimreduce1 = "scvis",
                                   dimreduce2 = "scvis",
                                   name1 = "follicular",
                                   name2 = "RLN",
                                   n_dims = 2,
                                   keep_cols = c("malignant_status_manual", "celltype_full"))



# Save outputs
saveRDS(scvis_map_RLN, args$out_rln)
saveRDS(scvis_combined, args$out_merged)

cat("Completed.\n")



