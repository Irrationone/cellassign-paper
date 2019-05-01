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
parser$add_argument('--hgsc_gene_list', metavar = 'FILE', type = 'character',
                    help="Path to HGSC gene list")
parser$add_argument('--hgsc_follicular_gene_list_combined', type = 'character',
                    help="Path to combined follicular + HGSC gene list")
parser$add_argument('--mcparland_gene_list_raw', metavar = 'FILE', type = 'character',
                    help="Path to McParland et al. raw gene list")
parser$add_argument('--mcparland_gene_list_revised', metavar = 'FILE', type = 'character',
                    help="Path to McParland et al. revised gene list")
parser$add_argument('--mcparland_gene_list_3types_revised', metavar = 'FILE', type = 'character',
                    help="Path to McParland et al. revised gene list for 3 types")
parser$add_argument('--liver_panglaodb_gene_list', metavar = 'FILE', type = 'character',
                    help="Path to PanglaoDB gene list")
parser$add_argument('--cellbench_gene_list_20', metavar = 'FILE', type = 'character',
                    help="Path to CellBench marker gene list (20)")
parser$add_argument('--cellbench_gene_list_30', metavar = 'FILE', type = 'character',
                    help="Path to CellBench marker gene list (30)")
parser$add_argument('--cellbench_gene_list_50', metavar = 'FILE', type = 'character',
                    help="Path to CellBench marker gene list (50)")
parser$add_argument('--cellbench_gene_list_73', metavar = 'FILE', type = 'character',
                    help="Path to CellBench marker gene list (73)")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()

koh_rho_path <- args$koh_rho
koh_rho <- read.table(koh_rho_path, sep = "\t", header = TRUE, row.names = 1) %>%
  tibble::rownames_to_column(var = "Gene")

follicular_gene_list_specific_path <- args$follicular_gene_list_specific
follicular_gene_list_lambda_kappa_path <- args$follicular_gene_list_lambda_kappa
hgsc_gene_list_path <- args$hgsc_gene_list
hgsc_follicular_gene_list_combined_path <- args$hgsc_follicular_gene_list_combined
mcparland_gene_list_raw_path <- args$mcparland_gene_list_raw
mcparland_gene_list_revised_path <- args$mcparland_gene_list_revised
mcparland_gene_list_3types_revised_path <- args$mcparland_gene_list_3types_revised
liver_panglaodb_gene_list_path <- args$liver_panglaodb_gene_list
cellbench_gene_list_20_path <- args$cellbench_gene_list_20
cellbench_gene_list_30_path <- args$cellbench_gene_list_30
cellbench_gene_list_50_path <- args$cellbench_gene_list_50
cellbench_gene_list_73_path <- args$cellbench_gene_list_73

gene_list_paths <- c(follicular_gene_list_specific_path,
                     follicular_gene_list_lambda_kappa_path,
                     hgsc_gene_list_path,
                     hgsc_follicular_gene_list_combined_path,
                     mcparland_gene_list_raw_path,
                     mcparland_gene_list_revised_path,
                     mcparland_gene_list_3types_revised_path,
                     liver_panglaodb_gene_list_path,
                     cellbench_gene_list_20_path,
                     cellbench_gene_list_30_path,
                     cellbench_gene_list_50_path,
                     cellbench_gene_list_73_path)
gene_list_names <- c("Koh et al.",
                     "FL celltype",
                     "FL light chain",
                     "HGSC celltype",
                     "HGSC and FL combined",
                     "McParland et al. liver",
                     "McParland et al. liver revised",
                     "McParland et al. revised 3 types",
                     "Liver PanglaoDB canonical",
                     "CellBench 20 genes",
                     "CellBench 30 genes",
                     "CellBench 50 genes",
                     "CellBench all 73 genes")

rho_dfs <- lapply(seq_along(gene_list_paths), function(i) {
  f <- gene_list_paths[i]
  marker_list <- read_yaml(f)
  rho <- marker_list_to_mat(marker_list, include_other = TRUE)
  rho_df <- as.data.frame(rho) %>%
    tibble::rownames_to_column(var = "Gene")
  return(rho_df)
})

rho_dfs <- c(list(as.data.frame(koh_rho)), rho_dfs)

# Write supplemental table xlsx
excel_export(rho_dfs,
             file = args$outfname,
             table_names = gene_list_names)

cat("Completed.\n")


