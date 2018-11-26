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
parser$add_argument('--deprob_methods', type='character', nargs ='+',
                    help="Clustering methods to use for DE prob analysis.")
parser$add_argument('--eval_measures', type='character', nargs ='+',
                    help="List of evaluation measures to show.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()

deprob_result_dir_naivecd8_naivecd4 <- args$deprob_result_dir_naivecd8_naivecd4
deprob_result_dir_b_cd8 <- args$deprob_result_dir_b_cd8
wrongmarker_result_dir_naivecd8_naivecd4 <- args$wrongmarker_result_dir_naivecd8_naivecd4
wrongmarker_result_dir_b_cd8 <- args$wrongmarker_result_dir_b_cd8
deprob_methods <- unlist(args$deprob_methods)
eval_measures <- unlist(args$eval_measures)

factor_orderings <- factor_orders()

naivecd8_naivecd4_measures <- load_annotation_files(deprob_result_dir_naivecd8_naivecd4, pattern = "*_eval_measures.tsv")
b_cd8_measures <- load_annotation_files(deprob_result_dir_b_cd8, pattern = "*_eval_measures.tsv")
naivecd8_naivecd4_deltas <- load_annotation_files(deprob_result_dir_naivecd8_naivecd4, pattern = "*_delta_compare.tsv")
b_cd8_deltas <- load_annotation_files(deprob_result_dir_b_cd8, pattern = "*_delta_compare.tsv")

# DE prob tables
naivecd8_naivecd4_measures <- naivecd8_naivecd4_measures %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "de") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "Naive CD8 vs. Naive CD4")

b_cd8_measures <- b_cd8_measures %>%
  dplyr::filter(clustering_method %in% deprob_methods,
                mapping_type == "de") %>%
  dplyr::select(c("seed", "gene_set", "num_groups", "num_cells",
                  "group_probs", "clustering_method",
                  "test_proportion", "sim_model", "de_prob",
                  eval_measures)) %>%
  dplyr::mutate(param_settings = "B vs. CD8")

eval_measures_merged <- plyr::rbind.fill(
  naivecd8_naivecd4_measures,
  b_cd8_measures
)

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

naivecd8_naivecd4_wm_eval_measures <- load_annotation_files(wrongmarker_result_dir_naivecd8_naivecd4, pattern = "*_eval_measures.tsv")
b_cd8_wm_eval_measures <- load_annotation_files(wrongmarker_result_dir_b_cd8, pattern = "*_eval_measures.tsv")

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
                  rvals_merged,
                  wm_eval_measures_merged),
             file = args$outfname,
             table_names = c("DE probability",
                             "Delta correspondence",
                             "Wrong marker proportion"))

cat("Completed.\n")


