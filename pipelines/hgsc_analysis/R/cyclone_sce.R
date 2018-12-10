#' Cell cycle prediction with cyclone on SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Run cyclone on SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--ncpus', type='double',
                    help="Number of cores to use.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cell cycle assignments.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

# Run cyclone
assignments <- cyclone(sce, hs.pairs, gene.names=rowData(sce)$ID, 
                       min.iter = 10, verbose = TRUE, BPPARAM = MulticoreParam(args$ncpus))

# Write outputs
saveRDS(assignments, file = args$outfname)

cat("Completed.\n")


