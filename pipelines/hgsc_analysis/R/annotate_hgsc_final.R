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
parser$add_argument('--unsupervised_epithelial', metavar='FILE', type='character',
                    help="Path to unsupervised epithelial assignments")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Read in assignments
cyclone_results <- readRDS(args$cyclone)
epithelial_assignments <- fread(args$unsupervised_epithelial)

# Add cell cycle information
sce@colData <- bind_cols(sce@colData %>% as.data.frame, 
                         cyclone_results$normalized.scores) %>% 
  DataFrame(check.names = FALSE)
sce <- sce %>%
  scater::mutate(Cell_Cycle = cyclone_results$phases)

# Rename cluster names
rename_cols <- function(df, prefix = '') {
  df <- df %>%
    dplyr::mutate_at(vars(contains('cluster')),
                     funs(factor(.))) %>%
    dplyr::rename_at(vars(contains('cluster')), 
                     funs(paste0(prefix, .)))
  return(df)
}
epithelial_renamed <- rename_cols(epithelial_assignments, prefix = "epithelial_")

# Add unsupervised clustering columns
sce@colData <- sce@colData %>% 
  data.frame(check.names = FALSE) %>%
  dplyr::left_join(epithelial_renamed) %>%
  DataFrame(check.names = FALSE)

# Save annotated SCE
saveRDS(sce, args$outfname)

cat("Completed.\n")



