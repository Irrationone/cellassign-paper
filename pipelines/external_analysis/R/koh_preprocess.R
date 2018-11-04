# Preprocess Koh data

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(Matrix)
library(DuoClustering2018)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Preprocess Koh data.")
parser$add_argument('--koh_celltypes', type = 'character', nargs='+',
                    help="Celltypes to use in Koh analysis (i.e. intersection of bulk and scRNAseq data)")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for Koh SCE")

koh_celltypes <- unlist(args$koh_celltypes)

# Use DuoClustering's version
sce_koh <- DuoClustering2018::sce_filteredExpr10_Koh()

sce_koh <- sce_koh %>%
  scater::mutate(celltype = str_replace_all(phenoid, "^H7(_derived_|_dreived_)?", "")) 

sce_koh_filtered <- sce_koh %>% 
  scater::filter(celltype %in% koh_celltypes)

# Normalization
sce_koh_normalized <- normalize(sce_koh_filtered)
sce_koh_normalized <- runPCA(sce_koh_normalized, ntop = 500, 
                             ncomponents = 50, exprs_values = "logcounts")
sce_koh_normalized <- runTSNE(sce_koh_normalized, use_dimred = "PCA", 
                              n_dimred = 50, ncomponents = 2)
sce_koh_normalized <- runUMAP(sce_koh_normalized, use_dimred = "PCA", 
                              n_dimred = 50, ncomponents = 2)

rownames(sce_koh_normalized) <- rownames(sce_koh_normalized)  %>%
  str_replace_all("\\..*", "")

# Save outputs
saveRDS(sce_koh_normalized, file = args$outfname)

cat("Completed.\n")

