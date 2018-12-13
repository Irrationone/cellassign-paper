#' Create SingleCellExperiment from UMI count files

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create SingleCellExperiment from UMI count files")
parser$add_argument('--sample_names', type='character', nargs = '+',
                    help="Sample name labels")
parser$add_argument('--umicount_files', metavar='FILE', type='character', nargs = '+',
                    help="UMI count files")
parser$add_argument('--patients', type='character', nargs = '+',
                    help="Patient labels")
parser$add_argument('--timepoints', type='character', nargs = '+',
                    help="Timepoint labels")
parser$add_argument('--sites', type='character', nargs = '+',
                    help="Anatomic sites")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for SingleCellExperiment.")
args <- parser$parse_args()

umicount_files <- unlist(args$umicount_files)

sces <- lapply(seq_along(umicount_files), function(i) {
  f <- umicount_files[i]
  sample_name <- args$sample_names[[i]]
  cnts <- fread(f) %>% 
    as.data.frame %>% 
    tibble::column_to_rownames(var = "CellId")
  cnts <- t(cnts)
  
  rowdat <- data.frame(Symbol=rownames(cnts), ID=rownames(cnts))
  rownames(rowdat) <- rowdat$Symbol
  
  coldat <- data.frame(Sample=sample_name, 
                       Barcode=colnames(cnts),
                       sample_barcode=paste0(sample_name, "_", colnames(cnts)),
                       dataset=sample_name,
                       timepoint=args$timepoints[[i]],
                       site=args$sites[[i]],
                       patient=args$patients[[i]])
  rownames(coldat) <- coldat$sample_barcode
  colnames(cnts) <- rownames(coldat)
  
  sce <- SingleCellExperiment(
    assays = list(counts = cnts),
    rowData = rowdat %>% DataFrame(),
    colData = coldat %>% DataFrame()
  )
  return(sce)
})

sce <- do.call('cbind', sces)


# Get ensembl gene IDs of ribosomal genes. Mitochondrial has already been filtered out
ribo_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)")]) %>%
  get_ensembl_id(sce)

# Calculate basic QC stats
sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
                            list(ribosomal=ribo_genes))

# Write outputs
saveRDS(sce, file = args$outfname)

cat("Completed.\n")


