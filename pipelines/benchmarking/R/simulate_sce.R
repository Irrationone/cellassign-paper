# Simulate SCE file

library(tidyverse)
library(scater)
library(scran)
library(R.utils)
library(methods)
library(argparse)

library(splatter)
library(cellassign.utils)


parser <- ArgumentParser(description = "Run CellAssign on a SingleCellExperiment")
parser$add_argument('--num_groups', type='integer',
                    help="Number of groups", default = 3)
parser$add_argument('--num_batches', type='integer',
                    help="Number of batches", default = 1)
parser$add_argument('--num_cells', type='integer',
                    help="Number of cells", default = 1000)
parser$add_argument('--num_genes', type='integer',
                    help="Number of genes", default = NULL)
parser$add_argument('--group_probs', type = 'character',
                    help="Group probabilities", default = "0.33,0.33,0.34")
parser$add_argument('--batch_probs', type = 'character',
                    help="Batch probabilities", default = "1")
parser$add_argument('--de_facscale', type = 'double',
                    help="DE variance", default = 0.02)
parser$add_argument('--de_facloc', type = 'double',
                    help="DE mean", default = 0.04)
parser$add_argument('--de_nu', type = 'double',
                    help="DE dispersion", default = 1.5)
parser$add_argument('--de_min', type = 'double',
                    help="Minimum logFC", default = 1)
parser$add_argument('--de_max', type = 'double',
                    help="Maximum logFC", default = 1000)
parser$add_argument('--batch_facscale', type = 'double',
                    help="Batch variance", default = 1)
parser$add_argument('--batch_facloc', type = 'double',
                    help="Batch mean", default = 1)
parser$add_argument('--de_prob', type = 'double',
                    help="DE probability", default = 0.25)
parser$add_argument('--down_prob', type = 'double',
                    help="Down regulation probability", default = 0.5)
parser$add_argument('--seed', type = 'integer',
                    help="Random seed", default = 1483)
parser$add_argument('--base_sce_param', type = 'character', metavar = 'FILE',
                    help="Base SCE param file.")
parser$add_argument('--sim_model', type = 'character',
                    help="Simulation model", default = "splat")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for simulated SCE")
args <- parser$parse_args()


num_groups <- args$num_groups
num_batches <- args$num_batches
num_cells <- args$num_cells
num_genes <- args$num_genes
group_probs <- as.numeric(strsplit(args$group_probs, ",")[[1]])
batch_probs <- as.numeric(strsplit(args$batch_probs, ",")[[1]])
de_facscale <- args$de_facscale
de_facloc <- args$de_facloc
de_nu <- args$de_nu
de_min <- args$de_min
de_max <- args$de_max
batch_facscale <- args$batch_facscale
batch_facloc <- args$batch_facloc
de_prob <- args$de_prob
down_prob <- args$down_prob
seed <- args$seed
base_sce_param_file <- args$base_sce_param
sim_model <- args$sim_model

param_df <- data.frame(num_groups=num_groups,
                       num_batches=num_batches,
                       num_cells=num_cells,
                       num_genes=num_genes,
                       group_probs=args$group_probs,
                       batch_probs=args$batch_probs,
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
                       base_sce_param=base_sce_param_file,
                       sim_model=sim_model, 
                       stringsAsFactors = FALSE)

model_params <- readRDS(base_sce_param_file)

group_probs <- group_probs/sum(group_probs)
batch_probs <- batch_probs/sum(batch_probs)

batch_cell_counts <- round(batch_probs * num_cells, 0)

num_cells <- sum(batch_cell_counts)
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

if (!is.null(num_genes)) {
  model_params_modified@nGenes <- num_genes
}

cat("Simulating ...\n")

sce_sim <- splatSimulate(params = model_params_modified, method = "groups", verbose = TRUE)

## Compute size factors
sce_sim <- calculateQCMetrics(sce_sim, exprs_values = "counts")
norm_factors <- edgeR::calcNormFactors(as.matrix(counts(sce_sim)), 
                                       method = "TMM")
lib_size_factors <- colSums(as.matrix(counts(sce_sim)))
sizeFactors(sce_sim) <- norm_factors * lib_size_factors/mean(norm_factors * 
                                                               lib_size_factors)

sce_sim@metadata$param_df <- param_df

saveRDS(sce_sim, args$outfname)

cat("Completed.\n")

