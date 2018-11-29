# Collate results of benchmarking

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(microbenchmark)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Collate results from benchmarking")
parser$add_argument('--timing_files', metavar='FILE', type='character', nargs ='+',
                    help="Timing output files")
parser$add_argument('--benchmark_config', metavar='FILE', type='character',
                    help="Benchmark config")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for stats.")
args <- parser$parse_args()

timing_files <- unlist(args$timing_files)

benchmark_config <- yaml::read_yaml(args$benchmark_config)
simulation_settings <- benchmark_config$simulation_settings
evaluation_settings <- benchmark_config$evaluation_settings

time_results <- plyr::rbind.fill(lapply(timing_files, function(f) {
  times <- readRDS(f)
  paths <- strsplit(f, "/")[[1]]
  params <- paths[(length(paths)-3):length(paths)]
  df <- data.frame(sim_seed=params[1],
                   eval_seed=tools::file_path_sans_ext(params[4]),
                   as.data.frame(simulation_settings[[params[2]]]),
                   as.data.frame(evaluation_settings[[params[3]]]),
                   time=times$time/1e9)
  return(df)
}))

# Write outputs

write.table(time_results, file = args$outfname,
            row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")

