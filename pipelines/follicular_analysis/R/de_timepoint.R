# Differential expression (grouped by timepoint)

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(gage)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Perform DE tests on an SCE for a particular celltype(s)")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--comparison', type='character',
                    help="Comparison type for DE")
parser$add_argument('--celltypes', type='character', nargs='+',
                    help="Celltypes to include")
parser$add_argument('--gene_set_file', type='character',
                    help="Gene set file to use")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

de_tables <- gage_analysis(sce = sce %>% scater::filter(celltype_full %in% unlist(args$celltypes)), 
                           comparison_type = "timepoint",
                           filter_ribo = filter_ribo,
                           filter_mito = filter_mito,
                           gene_set_path = args$gene_set_file, 
                           rowdat = rowData(sce))

# Save DE result (pathway and gene tables)
saveRDS(de_tables, args$outfname)

cat("Completed.\n")



