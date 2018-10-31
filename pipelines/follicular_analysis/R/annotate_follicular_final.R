# Annotate FL SCE with unsupervised cluster assignments and cyclone results

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
parser$add_argument('--unsupervised_t', metavar='FILE', type='character',
                    help="Path to unsupervised T assignments")
parser$add_argument('--unsupervised_malignant', metavar='FILE', type='character',
                    help="Path to unsupervised malignant B assignments")
parser$add_argument('--unsupervised_b', metavar='FILE', type='character',
                    help="Path to unsupervised nonmalignant B assignments")
parser$add_argument('--cyclone', metavar='FILE', type='character',
                    help="Path to cyclone results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated SCE.")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Read in assignments
t_assignments <- fread(args$unsupervised_t)
malignant_assignments <- fread(args$unsupervised_malignant)
b_assignments <- fread(args$unsupervised_b)
cyclone_results <- readRDS(args$cyclone)

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

t_renamed <- rename_cols(t_assignments, "t_")
malignant_renamed <- rename_cols(malignant_assignments, prefix = "malignant_")
b_renamed <- rename_cols(b_assignments, prefix = "b_")

# Add unsupervised clustering columns
sce@colData <- sce@colData %>% 
  data.frame(check.names = FALSE) %>%
  dplyr::left_join(t_renamed) %>%
  dplyr::left_join(malignant_renamed) %>%
  dplyr::left_join(b_renamed) %>% 
  DataFrame(check.names = FALSE)

# Save annotated SCE
saveRDS(sce, args$outfname)

cat("Completed.\n")



