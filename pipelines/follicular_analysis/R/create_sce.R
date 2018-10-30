#' Create SingleCellExperiment from filtered matrices

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create SingleCellExperiment from filtered matrices")
parser$add_argument('--sample_names', type='character', nargs = '+',
                    help="Sample name labels")
parser$add_argument('--filtered_matrices', metavar='DIR', type='character', nargs = '+',
                    help="Filtered matrix parent directories")
parser$add_argument('--patients', type='character', nargs = '+',
                    help="Patient labels")
parser$add_argument('--timepoints', type='character', nargs = '+',
                    help="Timepoint labels")
parser$add_argument('--progression', type='character', nargs = '+',
                    help="Progression type (primary, transformed, or progressed)")
parser$add_argument('--patient_progression', type='character', nargs = '+',
                    help="Overall progression status of patient (transformed or progressed)")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for SingleCellExperiment.")
args <- parser$parse_args()

filtered_matrix_dirs <- unlist(args$filtered_matrices)
names(filtered_matrix_dirs) <- unlist(args$sample_names)

sce <- DropletUtils::read10xCounts(unname(filtered_matrix_dirs))
sce <- sce %>% scater::mutate(
  dataset = str_extract(Sample, "FL[0-9]+[A-Z0-9_]+"),
  sample_barcode = paste(dataset, Barcode, sep = "_")
)

metadata_table <- data.frame(
  Sample=unname(filtered_matrix_dirs),
  dataset=unlist(args$sample_names),
  patient=unlist(args$patients),
  timepoint=unlist(args$timepoints),
  progression_status=unlist(args$progression),
  patient_progression=unlist(args$patient_progression),
  stringsAsFactors = FALSE
)

# Add metadata to sce
colData(sce) <- colData(sce) %>%
  as.data.frame %>%
  dplyr::left_join(metadata_table, by = 'Sample') %>%
  DataFrame()

# Get ensembl gene IDs of mitochondrial and ribosomal genes
mito_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^MT\\-")]) %>% 
  get_ensembl_id(sce)

ribo_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)")]) %>%
  get_ensembl_id(sce)

# Calculate basic QC stats
sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
                            list(mitochondrial=mito_genes, ribosomal=ribo_genes))

# Write outputs
saveRDS(sce, file = args$outfname)

cat("Completed.\n")


