
library(tidyverse)
library(scater)
library(scran)
library(R.utils)

library(splatter)
library(cellassign.utils)
library(methods)
library(argparse)

parser <- ArgumentParser(description = "Create simulated SCE.")
parser$add_argument('--params', metavar='FILE', type='character',
                    help="Splatter base params object")
parser$add_argument('--outfname', metavar='FILE', type = 'character',
                    help="Output filename for SCE.")
parser$add_argument('--seed', type = 'integer',
                    help="Random seed")
parser$add_argument('--sim_model', type = 'character',
                    help="Simulation model", default = 'splat')
parser$add_argument('--num_cells', type = 'integer',
                    help="Number of cells", default = 1000)
parser$add_argument('--num_groups', type = 'integer',
                    help="Number of groups", default = 2)
parser$add_argument('--num_batches', type = 'integer',
                    help="Number of batches")
parser$add_argument('--group_probs', type = 'double', nargs = '+',
                    help="Group probabilities")
parser$add_argument('--batch_probs', type = 'double', nargs = '+',
                    help="Batch probabilities")
parser$add_argument('--de_facloc', type = 'double',
                    help="DE location parameter")
parser$add_argument('--de_facscale', type = 'double',
                    help="DE scale parameter")
parser$add_argument('--de_nu', type = 'double',
                    help="DE nu parameter")
parser$add_argument('--de_min', type = 'double',
                    help="Min FC")
parser$add_argument('--de_max', type = 'double',
                    help="Max FC")
parser$add_argument('--batch_facloc', type = 'double',
                    help="Batch location parameter")
parser$add_argument('--batch_facscale', type = 'double',
                    help="Batch scale parameter")
parser$add_argument('--de_prob', type = 'double',
                    help="DE probability")
parser$add_argument('--down_prob', type = 'double',
                    help="Downregulation probability")
args <- parser$parse_args()

group_probs <- unlist(args$group_probs)
batch_probs <- unlist(args$batch_probs)
num_groups <- args$num_groups
num_batches <- args$num_batches
num_cells <- args$num_cells
de_facloc <- args$de_facloc
de_facscale <- args$de_facscale
de_nu <- args$de_nu
de_min <- args$de_min
de_max <- args$de_max
batch_facloc <- args$batch_facloc
batch_facscale <- args$batch_facscale
de_prob <- args$de_prob
down_prob <- args$down_prob
seed <- args$seed

stopifnot(length(group_probs) == num_groups)
stopifnot(length(batch_probs) == num_batches)

model_params <- readRDS(args$params)

# Create parameter dataframe
param_df <- data.frame(num_groups=num_groups,
                       num_batches=num_batches,
                       num_cells=num_cells,
                       group_probs=paste(group_probs, collapse = ","),
                       batch_probs=paste(batch_probs, collapse = ","),
                       de_facscale=de_facscale,
                       de_facloc=de_facloc,
                       de_nu=de_nu,
                       de_min=de_min,
                       de_max=de_max,
                       batch_facscale=batch_facscale,
                       batch_facloc=batch_facloc,
                       de_prob=de_prob,
                       down_prob=down_prob,
                       seed=seed,
                       base_sce_param=args$params,
                       sim_model=args$sim_model, stringsAsFactors = FALSE)

# Renormalize probabilities
group_probs <- group_probs/sum(group_probs)
batch_probs <- batch_probs/sum(batch_probs)

batch_cell_counts <- round(batch_probs * num_cells, 0)
num_cells <- sum(batch_cell_counts)

# Overwrite model params
model_params_modified <- model_params
model_params_modified@nCells <- num_cells
model_params_modified@de.prob <- de_prob
model_params_modified@de.downProb <- down_prob
model_params_modified@nBatches <- num_batches
model_params_modified@nGroups <- num_groups
model_params_modified@batchCells <- batch_cell_counts
model_params_modified@group.prob <- group_probs
model_params_modified@de.facLoc <- de_facloc
model_params_modified@de.facScale <- de_facscale
model_params_modified@de.facNu <- de_nu
model_params_modified@de.max <- de_max
model_params_modified@de.min <- de_min
model_params_modified@batch.facLoc <- batch_facloc
model_params_modified@batch.facScale <- batch_facscale
model_params_modified@seed <- seed

sce_sim <- splatSimulate(params = model_params_modified, method = "groups", verbose = TRUE)

sce_sim <- normalize_sce_sim(sce_sim, run_dimred = TRUE)

sce_sim@metadata$param_df <- param_df

# Write out SCE
saveRDS(sce_sim, args$outfname)

