# Annotate HGSC SCE with cyclone results

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Annotate SCE with unsupervised clustering and cyclone results")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--cyclone', metavar='FILE', type='character',
                    help="Path to cyclone results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Read in assignments
cyclone_results <- readRDS(args$cyclone)

# Add cell cycle information
sce@colData <- bind_cols(sce@colData %>% as.data.frame, 
                         cyclone_results$normalized.scores) %>% 
  DataFrame(check.names = FALSE)
sce <- sce %>%
  scater::mutate(Cell_Cycle = cyclone_results$phases)


# Save annotated SCE
saveRDS(sce, args$outfname)

cat("Completed.\n")



