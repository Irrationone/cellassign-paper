# Filter follicular SCE by a particular celltype and malignant status
# Necessary to ensure consistency between results derived from the resulting SCE

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Extract cells from follicular SCE")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--celltype', type = 'character',
                    help="Celltype to filter", 
                    choices = c("all", "T cell", "B cell"))
parser$add_argument('--celltype_full', type = 'character', nargs ="+",
                    help="Full celltypes to consider", default = NULL)
parser$add_argument('--celltype_probability', type='double',
                    help="Broad probability threshold to consider a cell a T cell", default = 0.9)
parser$add_argument('--malignancy_filter', type='character',
                    help="Filter for malignant/nonmalignant cells", 
                    choices = c("malignant", "nonmalignant", "all"),
                    default = "all")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for SCE")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
tcell_labels <- unlist(args$tcell_labels)

categorical_palettes <- cat_palettes()

# Filter for celltype
if (args$celltype != "all") {
  sce <- follicular_filter(sce, filter = args$celltype, probability = args$celltype_probability)
}

# Filter for full celltype 
if (!is.null(args$celltype_full)) {
  sce <- sce %>%
    scater::filter(celltype_full %in% unlist(args$celltype_full))
}

# Filter for malignant status
if (args$malignancy_filter != "all") {
  sce <- sce %>%
    scater::filter(malignant_status_manual == args$malignancy_filter)
}


# Dimensionality reduction (again)
sce <- runPCA(sce, ntop = 500, ncomponents = 50)
sce <- runTSNE(sce, ncomponents = 2, use_dimred = "PCA", n_dimred = 50)
sce <- runUMAP(sce, ncomponents = 2, use_dimred = "PCA", n_dimred = 50)

# Save outputs
saveRDS(sce, args$outfname)

cat("Completed.\n")




