# Supplemental table containing marker gene matrices

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
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create rho supplemental table")
parser$add_argument('--koh_rho', metavar='FILE', type='character',
                    help="Path to Koh rho")
parser$add_argument('--follicular_gene_list_specific', metavar = 'FILE', type = 'character',
                    help="Path to follicular gene list (specific)")
parser$add_argument('--follicular_gene_list_lambda_kappa', metavar = 'FILE', type = 'character',
                    help="Path to follicular gene list (lambda/kappa)")
parser$add_argument('--follicular_specific_other', type = 'character',
                    help="Include an other column in the marker gene matrix if using a marker_list")
parser$add_argument('--follicular_lk_other', type = 'character',
                    help="Include an other column in the marker gene matrix if using a marker_list")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()

koh_rho_path <- args$koh_rho
koh_rho <- read.table(koh_rho_path, sep = "\t", header = TRUE, row.names = 1)

follicular_gene_list_specific_path <- args$follicular_gene_list_specific
follicular_gene_list_lambda_kappa_path <- args$follicular_gene_list_lambda_kappa
follicular_specific_other <- as.logical(args$follicular_specific_other)
follicular_lk_other <- as.logical(args$follicular_lk_other)

gene_list_paths <- c(follicular_gene_list_specific_path,
                     follicular_gene_list_lambda_kappa_path)
include_other_settings <- c(follicular_specific_other,
                            follicular_lk_other)
gene_list_names <- c("Koh et al.",
                     "FL celltype",
                     "FL light chain")

rho_dfs <- lapply(seq_along(gene_list_paths), function(i) {
  f <- gene_list_paths[i]
  marker_list <- read_yaml(f)
  rho <- marker_list_to_mat(marker_list, include_other = include_other_settings[i])
  rho_df <- as.data.frame(rho)
  return(rho_df)
})

rho_dfs <- c(list(as.data.frame(koh_rho)), rho_dfs)

# Write supplemental table xlsx
excel_export(rho_dfs,
             file = args$outfname,
             table_names = gene_list_names)

cat("Completed.\n")


