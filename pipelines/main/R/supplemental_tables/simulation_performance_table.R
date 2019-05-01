# Supplemental table (simulation)

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(Matrix)
library(xlsx)
library(ImportExport)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create simulation supp table.")
parser$add_argument('--deprob_result_dir_naivecd8_naivecd4', metavar='DIR', type='character',
                    help="Path to simulation result directory for naive CD8+ vs naive CD4+")
parser$add_argument('--deprob_result_dir_b_cd8', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for B vs. CD8")
parser$add_argument('--wrongmarker_result_dir_naivecd8_naivecd4', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker (naive CD8+ vs naive CD4+")
parser$add_argument('--wrongmarker_result_dir_b_cd8', metavar = 'DIR', type = 'character',
                    help="Path to simulation result directory for wrong marker (B vs. CD8+)")
parser$add_argument('--novel_extra_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for novel/extra celltype analysis")
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--eval_measures', type='character', nargs ='+',
                    help="List of evaluation measures to show.")
parser$add_argument('--method_description', type='character', metavar='FILE',
                    help="Method description")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()

deprob_result_dir_naivecd8_naivecd4 <- args$deprob_result_dir_naivecd8_naivecd4
deprob_result_dir_b_cd8 <- args$deprob_result_dir_b_cd8
wrongmarker_result_dir_naivecd8_naivecd4 <- args$wrongmarker_result_dir_naivecd8_naivecd4
wrongmarker_result_dir_b_cd8 <- args$wrongmarker_result_dir_b_cd8
deprob_methods <- unlist(args$deprob_methods)
eval_measures <- unlist(args$eval_measures)
novel_extra_result_dir <- args$novel_extra_result_dir

method_description <- read_yaml(args$method_description)
method_metadata <- plyr::rbind.fill(lapply(names(method_description), function(i) {
  data.frame(method_type=i, clustering_method=method_description[[i]])
}))

factor_orderings <- factor_orders()

load_sim_performance <- function(deprob_result_dir, method_metadata) {
  deprob_result_files <- Sys.glob(file.path(deprob_result_dir, "evaluate", "*", "*", "*.tsv"))
  de_eval_measures <- plyr::rbind.fill(lapply(deprob_result_files, function(f) {
    df <- fread(f)
    if ("gene_set" %in% colnames(df)) {
      return(df)
    } else {
      feature_type <- str_extract(f, "(markers|full)")
      return(data.frame(fread(f), gene_set=feature_type))
    }
  }))
  de_eval_measures <- de_eval_measures %>%
    dplyr::left_join(method_metadata)
  return(de_eval_measures)
}

load_deltas <- function(deprob_result_dir) {
  delta_files <- Sys.glob(file.path(deprob_result_dir, "assign_celltypes_sce", "deltas", "*", "cellassign*.tsv"))
  de_deltas <- plyr::rbind.fill(lapply(delta_files, function(f) {
    fread(f)
  }))
  return(de_deltas)
}

naivecd8_naivecd4_measures_raw <- load_sim_performance(deprob_result_dir_naivecd8_naivecd4, method_metadata)
b_cd8_measures_raw <- load_sim_performance(deprob_result_dir_b_cd8, method_metadata)
naivecd8_naivecd4_deltas <- load_deltas(deprob_result_dir_naivecd8_naivecd4)
b_cd8_deltas <- load_deltas(deprob_result_dir_b_cd8)
#naivecd8_naivecd4_measures <- load_annotation_files(deprob_result_dir_naivecd8_naivecd4, pattern = "*_eval_measures.tsv")
#b_cd8_measures <- load_annotation_files(deprob_result_dir_b_cd8, pattern = "*_eval_measures.tsv")
#naivecd8_naivecd4_deltas <- load_annotation_files(deprob_result_dir_naivecd8_naivecd4, pattern = "*_delta_compare.tsv")
#b_cd8_deltas <- load_annotation_files(deprob_result_dir_b_cd8, pattern = "*_delta_compare.tsv")

novel_extra_measures <- load_sim_performance(novel_extra_result_dir, method_metadata)

# DE prob tables
naivecd8_naivecd4_measures <- naivecd8_naivecd4_measures_raw %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "de") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "Naive CD8 vs. Naive CD4")

b_cd8_measures <- b_cd8_measures_raw %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "de") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "B vs. CD8")

naivecd8_naivecd4_measures_correlation <- naivecd8_naivecd4_measures_raw %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "correlation") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "Naive CD8 vs. Naive CD4 (correlation)")

b_cd8_measures_correlation <- b_cd8_measures_raw %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "correlation") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "B vs. CD8 (correlation)")

eval_measures_merged <- plyr::rbind.fill(
  naivecd8_naivecd4_measures,
  b_cd8_measures,
  naivecd8_naivecd4_measures_correlation,
  b_cd8_measures_correlation
)

# Novel/extra celltype table

novel_extra_measures_filtered <- novel_extra_measures %>%
  dplyr::select(c("seed", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "sim_model", "de_prob", "n_data_types", "n_marker_types",
                  "data_subset_types", "marker_subset_types",
                  eval_measures))

# Delta tables (report correlation)

naivecd8_naivecd4_rvals <- compute_pvals_subsets(naivecd8_naivecd4_deltas,
                                                 facet_vars = c("de_prob", "clustering_method"),
                                                 formula = ~ true_delta + inferred_delta,
                                                 corfun = cor.test,
                                                 output = "estimate") %>%
  dplyr::mutate(param_settings = "Naive CD8 vs. Naive CD4")

b_cd8_rvals <- compute_pvals_subsets(b_cd8_deltas,
                                     facet_vars = c("de_prob", "clustering_method"),
                                     formula = ~ true_delta + inferred_delta,
                                     corfun = cor.test,
                                     output = "estimate") %>%
  dplyr::mutate(param_settings = "B vs. CD8")

rvals_merged <- plyr::rbind.fill(naivecd8_naivecd4_rvals,
                                 b_cd8_rvals) %>%
  dplyr::rename(correlation_coefficient=estimate)

# Wrong marker tables

naivecd8_naivecd4_wm_eval_measures <- load_sim_performance(wrongmarker_result_dir_naivecd8_naivecd4, method_metadata)
b_cd8_wm_eval_measures <- load_sim_performance(wrongmarker_result_dir_b_cd8, method_metadata)
#naivecd8_naivecd4_wm_eval_measures <- load_annotation_files(wrongmarker_result_dir_naivecd8_naivecd4, pattern = "*_eval_measures.tsv")
#b_cd8_wm_eval_measures <- load_annotation_files(wrongmarker_result_dir_b_cd8, pattern = "*_eval_measures.tsv")

naivecd8_naivecd4_wm_eval_measures <- naivecd8_naivecd4_wm_eval_measures %>% 
  dplyr::select(c("clustering_method", "max_genes", "wrong_marker_proportion", eval_measures)) %>%
  dplyr::mutate(param_settings = "Naive CD8 vs. Naive CD4")

b_cd8_wm_eval_measures <- b_cd8_wm_eval_measures %>% 
  dplyr::select(c("clustering_method", "max_genes", "wrong_marker_proportion", eval_measures)) %>%
  dplyr::mutate(param_settings = "B vs. CD8")

wm_eval_measures_merged <- plyr::rbind.fill(naivecd8_naivecd4_wm_eval_measures,
                                            b_cd8_wm_eval_measures)

# Write supplemental table xlsx
excel_export(list(eval_measures_merged,
                  novel_extra_measures_filtered,
                  rvals_merged,
                  wm_eval_measures_merged),
             file = args$outfname,
             table_names = c("DE probability",
                             "Missing and extra cell types",
                             "Delta correspondence",
                             "Wrong marker proportion"))

cat("Completed.\n")


