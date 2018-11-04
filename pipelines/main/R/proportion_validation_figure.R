# Figure for proportion validation

# Duo validation
# HGSC validation

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(pheatmap)
library(Matrix)
library(DuoClustering2018)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create proportion validation figure.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

# Koh analysis

common_celltypes <- c("hESC", "APS", "MPS", "DLL1pPXM", "ESMT", "Sclrtm", "D5CntrlDrmmtm",
                      "D2LtM") #"D3GARPpCrdcM")
koh_filtered <- koh %>% 
  scater::filter(celltype %in% common_celltypes)

get_ensembl_id <- function(x, rowdat, symbol_name = "symbol", gene_name = "gene") {
  rowdat <- as.data.frame(rowdat)
  df_as_map(rowdat %>% subset(x %in% rowdat[,symbol_name]), x, from = symbol_name, to = gene_name)
}

mito_genes <- as.character(rowData(koh_filtered)$symbol[str_detect(rowData(koh_filtered)$symbol, "^MT\\-")]) %>% 
  get_ensembl_id(rowData(koh_filtered))

ribo_genes <- as.character(rowData(koh_filtered)$symbol[str_detect(rowData(koh_filtered)$symbol, "^RP(L|S)")]) %>%
  get_ensembl_id(rowData(koh_filtered))

koh_filtered <- calculateQCMetrics(koh_filtered, exprs_values = "counts", feature_controls =
                                     list(mitochondrial=mito_genes, ribosomal=ribo_genes))

koh_filtered <- filter_cells(koh_filtered, nmads = Inf, type = "lower", log = TRUE, max_mito = 30, max_ribo = 50)

norm_factors <- edgeR::calcNormFactors(as.matrix(counts(koh_filtered)), 
                                       method = "TMM")
lib_size_factors <- colSums(as.matrix(counts(koh_filtered)))
sizeFactors(koh_filtered) <- norm_factors * lib_size_factors/mean(norm_factors * 
                                                                    lib_size_factors)
koh_normalized <- normalize(koh_filtered)
koh_normalized <- runPCA(koh_normalized, ntop = 500, 
                         ncomponents = 50, exprs_values = "logcounts")
koh_normalized <- runTSNE(koh_normalized, use_dimred = "PCA", 
                          n_dimred = 50, ncomponents = 2)



