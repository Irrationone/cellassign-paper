# Create a list of all the stats for the paper

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(ggrepel)
library(Matrix)
library(rjson)
library(yaml)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create compiled stats file.")
parser$add_argument('--sce_follicular', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for FL")
parser$add_argument('--sce_follicular_bcells', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for FL")
parser$add_argument('--sce_follicular_tcells', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for FL")
parser$add_argument('--sce_follicular_raw', metavar='FILE', type='character',
                    help="Path to raw SingleCellExperiment RDS")
parser$add_argument('--follicular_tcell_labels', nargs='+',
                    type='character', help='Follicular T cell types')
parser$add_argument('--patient_progression', metavar='FILE', type='character',
                    help="Patient events")
parser$add_argument('--conda_env',
                    type='character', help='Conda environment', default = 'r-tensorflow')
parser$add_argument('--deprob_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for DE prob")
parser$add_argument('--delta_deprobs', type='double', nargs ='+',
                    help="DE probs to show delta plots for.")
parser$add_argument('--method_description', metavar='FILE', type='character',
                    help="Description of competing methods")
parser$add_argument('--sce_follicular_nonmalignant_bcell', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--cellassign_lambda_kappa', metavar='FILE', type='character',
                    help="Path to CellAssign lambda/kappa results")
parser$add_argument('--de_timepoint_dir', type='character',
                    help="DE results for timepoint comparisons")
parser$add_argument('--de_timepoint_fgsea_dir', type='character',
                    help="DE results for timepoint comparisons (using FGSEA and hallmark gene set)")
parser$add_argument('--de_malignant_timepoint_dir', type='character',
                    help="DE results for malignant-timepoint interaction comparisons")
parser$add_argument('--koh_annotated', metavar='FILE', type='character',
                    help="Annotated SCE of Koh")
parser$add_argument('--koh_rho', type='character', metavar='FILE',
                    help="Koh et al. rho matrix")
parser$add_argument('--cellassign_koh', metavar='FILE', type='character',
                    help="CellAssign fit to Koh")
parser$add_argument('--scina_koh', metavar='FILE', type='character',
                    help="SCINA fit to Koh")
parser$add_argument('--sce_hgsc', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for HGSC")
parser$add_argument('--sce_hgsc_raw', metavar='FILE', type='character',
                    help="Path to raw SingleCellExperiment RDS for HGSC")
parser$add_argument('--de_site', metavar='DIR', type='character',
                    help="Directory to DE site results")
parser$add_argument('--de_site_fgsea', metavar='DIR', type='character',
                    help="Directory to DE site fgsea results")
parser$add_argument('--de_epithelial', metavar='DIR', type='character',
                    help="Directory to DE epithelial clusters results")
parser$add_argument('--sce_liver', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for liver")
parser$add_argument('--cellassign_liver_3markers_3types', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_4markers_3types', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_allmarkers_3types', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_3markers_4types', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_3markers_alltypes_raw', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_3markers_alltypes', metavar='FILE', type='character')
parser$add_argument('--scina_liver_3markers_3types', metavar='FILE', type='character')
parser$add_argument('--scina_liver_4markers_3types', metavar='FILE', type='character')
parser$add_argument('--scina_liver_allmarkers_3types', metavar='FILE', type='character')
parser$add_argument('--scina_liver_allmarkers_3types_point1', metavar='FILE', type='character')
parser$add_argument('--scina_liver_3markers_4types', metavar='FILE', type='character')
parser$add_argument('--scina_liver_3markers_alltypes_raw', metavar='FILE', type='character')
parser$add_argument('--scina_liver_3markers_alltypes', metavar='FILE', type='character')
parser$add_argument('--cellassign_liver_panglaodb', metavar='FILE', type='character')
parser$add_argument('--liver_3_celltypes', type='character', nargs = '+')
parser$add_argument('--liver_4_celltypes', type='character', nargs = '+')
parser$add_argument('--liver_all_celltypes', type='character', nargs = '+')
parser$add_argument('--liver_panglaodb_markers', metavar='FILE', type='character',
                    help="PanglaoDB liver markers")
parser$add_argument('--sce_pure_merged', metavar='FILE', type='character')
parser$add_argument('--sce_pure_10x', metavar='FILE', type='character')
parser$add_argument('--sce_pure_celseq', metavar='FILE', type='character')
parser$add_argument('--sce_pure_dropseq', metavar='FILE', type='character')
parser$add_argument('--sce_mix_merged', metavar='FILE', type='character')
parser$add_argument('--cellassign_tian_pure_merged', metavar='FILE', type='character')
parser$add_argument('--cellassign_tian_pure_10x', metavar='FILE', type='character')
parser$add_argument('--cellassign_tian_pure_celseq', metavar='FILE', type='character')
parser$add_argument('--cellassign_tian_pure_dropseq', metavar='FILE', type='character')
parser$add_argument('--cellassign_tian_mix', metavar='FILE', type='character')
parser$add_argument('--sce_fl_hgsc_merged', metavar='FILE', type='character')
parser$add_argument('--cellassign_fl_hgsc', metavar='FILE', type='character')
parser$add_argument('--cellassign_fl_combination', metavar='FILE', type='character')
parser$add_argument('--cellassign_hgsc_combination', metavar='FILE', type='character')
parser$add_argument('--cellassign_fl_noother', metavar='FILE', type='character')
parser$add_argument('--cellassign_fl_top', metavar='FILE', type='character')
parser$add_argument('--cellassign_fl_bottom', metavar='FILE', type='character')
parser$add_argument('--cellassign_koh_noother', metavar='FILE', type='character')
parser$add_argument('--cellassign_koh_top', metavar='FILE', type='character')
parser$add_argument('--cellassign_koh_bottom', metavar='FILE', type='character')
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for stats file")
args <- parser$parse_args()

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

sce_follicular <- readRDS(args$sce_follicular)
sce_follicular_raw <- readRDS(args$sce_follicular_raw)
sce_follicular_nonmalignant_b <- readRDS(args$sce_follicular_nonmalignant_bcell)
sce_follicular_b <- readRDS(args$sce_follicular_bcells)
deprob_result_dir <- args$deprob_result_dir
delta_deprobs <- unlist(args$delta_deprobs)
cellassign_lambda_kappa_results <- readRDS(args$cellassign_lambda_kappa)

method_description <- read_yaml(args$method_description)
method_metadata <- plyr::rbind.fill(lapply(names(method_description), function(i) {
  data.frame(method_type=i, clustering_method=method_description[[i]])
}))

sce_follicular_tcell <- readRDS(args$sce_follicular_tcell)
follicular_tcell_labels <- unlist(args$follicular_tcell_labels)
sce_follicular_tcell <- sce_follicular_tcell %>%
  scater::filter(celltype_full %in% c(follicular_tcell_labels, "other")) %>%
  scater::mutate(celltype_full=factor(plyr::mapvalues(celltype_full,
                                                      "other",
                                                      "Unassigned")))
sce_follicular_cytotoxic <- sce_follicular_tcell %>%
  scater::filter(celltype_full %in% c("Cytotoxic T cells"))

de_timepoint_dir <- args$de_timepoint_dir
de_timepoint_fgsea_dir <- args$de_timepoint_fgsea_dir
de_malignant_timepoint_dir <- args$de_malignant_timepoint_dir

de_timepoint_files <- Sys.glob(file.path(de_timepoint_dir, "*", "*"))
de_timepoint_fgsea_files <- Sys.glob(file.path(de_timepoint_fgsea_dir, "malignant", "*"))
de_malignant_timepoint_files <- Sys.glob(file.path(de_malignant_timepoint_dir, "*", "*"))

koh_rho <- read.table(args$koh_rho, sep = "\t", header = TRUE, row.names = 1) %>%
  tibble::rownames_to_column(var = "Gene")
koh_annotated_path <- args$koh_annotated
cellassign_fit_koh <- readRDS(args$cellassign_koh)
scina_fit_koh <- readRDS(args$scina_koh)

# HGSC data

sce_hgsc <- readRDS(args$sce_hgsc)
sce_hgsc_raw <- readRDS(args$sce_hgsc_raw)
de_site_dir <- args$de_site
de_site_fgsea_dir <- args$de_site_fgsea
de_epithelial_dir <- args$de_epithelial

# Liver data

sce_liver <- readRDS(args$sce_liver)
liver_fits <- list(
  'cellassign_liver_threemarkers_threetypes'=readRDS(args$cellassign_liver_3markers_3types),
  'cellassign_liver_fourmarkers_threetypes'=readRDS(args$cellassign_liver_4markers_3types),
  'cellassign_liver_allmarkers_threetypes'=readRDS(args$cellassign_liver_allmarkers_3types),
  'cellassign_liver_threemarkers_fourtypes'=readRDS(args$cellassign_liver_3markers_4types),
  'cellassign_liver_threemarkers_alltypes_raw'=readRDS(args$cellassign_liver_3markers_alltypes_raw),
  'cellassign_liver_threemarkers_alltypes'=readRDS(args$cellassign_liver_3markers_alltypes),
  'scina_liver_threemarkers_threetypes'=readRDS(args$scina_liver_3markers_3types),
  'scina_liver_fourmarkers_threetypes'=readRDS(args$scina_liver_4markers_3types),
  'scina_liver_allmarkers_threetypes'=readRDS(args$scina_liver_allmarkers_3types),
  'scina_liver_allmarkers_threetypes_pointone'=readRDS(args$scina_liver_allmarkers_3types_point1),
  'scina_liver_threemarkers_fourtypes'=readRDS(args$scina_liver_3markers_4types),
  'scina_liver_threemarkers_alltypes_raw'=readRDS(args$scina_liver_3markers_alltypes_raw),
  'scina_liver_threemarkers_alltypes'=readRDS(args$scina_liver_3markers_alltypes),
  'cellassign_liver_panglaodb'=readRDS(args$cellassign_liver_panglaodb)
)

liver_data_celltypes <- list(
  'cellassign_liver_threemarkers_threetypes'=unlist(args$liver_3_celltypes),
  'cellassign_liver_fourmarkers_threetypes'=unlist(args$liver_3_celltypes),
  'cellassign_liver_allmarkers_threetypes'=unlist(args$liver_3_celltypes),
  'cellassign_liver_threemarkers_fourtypes'=unlist(args$liver_4_celltypes),
  'cellassign_liver_threemarkers_alltypes_raw'=unlist(args$liver_all_celltypes),
  'cellassign_liver_threemarkers_alltypes'=unlist(args$liver_all_celltypes),
  'scina_liver_threemarkers_threetypes'=unlist(args$liver_3_celltypes),
  'scina_liver_fourmarkers_threetypes'=unlist(args$liver_3_celltypes),
  'scina_liver_allmarkers_threetypes'=unlist(args$liver_3_celltypes),
  'scina_liver_allmarkers_threetypes_pointone'=unlist(args$liver_3_celltypes),
  'scina_liver_threemarkers_fourtypes'=unlist(args$liver_4_celltypes),
  'scina_liver_threemarkers_alltypes_raw'=unlist(args$liver_all_celltypes),
  'scina_liver_threemarkers_alltypes'=unlist(args$liver_all_celltypes),
  'cellassign_liver_panglaodb'=unlist(args$liver_all_celltypes)
)

liver_panglaodb_markers <- read_yaml(args$liver_panglaodb_markers)

# CellBench data

sce_pure_merged <- readRDS(args$sce_pure_merged)
sce_pure_10x <- readRDS(args$sce_pure_10x)
sce_pure_celseq <- readRDS(args$sce_pure_celseq)
sce_pure_dropseq <- readRDS(args$sce_pure_dropseq)
sce_mix_merged <- readRDS(args$sce_mix_merged)

cellbench_fits <- list(
  'cellassign_tian_pure_merged'=readRDS(args$cellassign_tian_pure_merged),
  'cellassign_tian_mix'=readRDS(args$cellassign_tian_mix),
  'cellassign_tian_pure_tenx'=readRDS(args$cellassign_tian_pure_10x),
  'cellassign_tian_pure_celseq'=readRDS(args$cellassign_tian_pure_celseq),
  'cellassign_tian_pure_dropseq'=readRDS(args$cellassign_tian_pure_dropseq)
)

cellbench_sets <- list(
  'cellassign_tian_pure_merged'=sce_pure_merged,
  'cellassign_tian_mix'=sce_mix_merged,
  'cellassign_tian_pure_tenx'=sce_pure_10x,
  'cellassign_tian_pure_celseq'=sce_pure_celseq,
  'cellassign_tian_pure_dropseq'=sce_pure_dropseq
)

# Combination/hierarchical data

sce_fl_hgsc_merged <- readRDS(args$sce_fl_hgsc_merged)
combination_fits <- list(
  'cellassign_fl_hgsc'=readRDS(args$cellassign_fl_hgsc),
  'cellassign_fl_combination'=readRDS(args$cellassign_fl_combination),
  'cellassign_hgsc_combination'=readRDS(args$cellassign_hgsc_combination),
  'cellassign_fl_noother'=readRDS(args$cellassign_fl_noother),
  'cellassign_fl_top'=readRDS(args$cellassign_fl_top),
  'cellassign_fl_bottom'=readRDS(args$cellassign_fl_bottom),
  'cellassign_koh_noother'=readRDS(args$cellassign_koh_noother),
  'cellassign_koh_top'=readRDS(args$cellassign_koh_top),
  'cellassign_koh_bottom'=readRDS(args$cellassign_koh_bottom)
)

## Delta stats

message("Performing DE prob analysis ...")

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

delta_files <- Sys.glob(file.path(deprob_result_dir, "assign_celltypes_sce", "deltas", "*", "cellassign*.tsv"))
de_deltas <- plyr::rbind.fill(lapply(delta_files, function(f) {
  fread(f)
}))

delta_table <- de_deltas %>% 
  dplyr::filter(de_prob %in% delta_deprobs) 

rvals <- compute_pvals_subsets(delta_table,
                               facet_vars = c("de_prob", "clustering_method"),
                               formula = ~ true_delta + inferred_delta,
                               corfun = cor.test,
                               output = "estimate")

min_delta_rval <- round(min(rvals$estimate), 3)

## Koh rho

message("Koh rho stats ...")

num_koh_markers <- nrow(koh_rho)

## FL cell counts

message("Follicular lymphoma stats ...")

fl_total_cell_count <- ncol(sce_follicular_raw)
fl_cell_counts <- as.list(table(sce_follicular_raw$patient))
fl_cell_counts_by_sample <- as.list(table(sce_follicular_raw$dataset))

## Lambda/kappa stats

message("Lambda/kappa stats ...")

sce_follicular_nonmalignant_b$class <- cellassign_lambda_kappa_results$cell_type
lk_nonother <- sum(sce_follicular_nonmalignant_b$class != "other")
lk_kappa <- sum(sce_follicular_nonmalignant_b$class == "IGKC")
lk_lambda <- sum(sce_follicular_nonmalignant_b$class == "IGLC")
lk_kappa_pct <- round(lk_kappa/lk_nonother * 100, 1)
lk_kappa_patient <- table((sce_follicular_nonmalignant_b %>% scater::filter(class == "IGKC"))$patient)
lk_lambda_patient <- table((sce_follicular_nonmalignant_b %>% scater::filter(class == "IGLC"))$patient)
lk_nonother_patient <- table((sce_follicular_nonmalignant_b %>% scater::filter(class != "other"))$patient)
lk_kappa_pct_patient <- round(lk_kappa_patient/lk_nonother_patient * 100, 1)

## Malignant marker pvals

message("Malignant marker stats ...")

malignant_markers <- c("BCL2", "BCL6")
malignant_marker_pvals <- lapply(de_malignant_timepoint_files, function(f) {
  de_res <- readRDS(f)
  pvals <- (de_res$malignant$gene %>% 
              dplyr::filter(Symbol %in% malignant_markers))$adj.P.Val
  return(pvals)
})
max_malignant_pval <- signif(max(unlist(malignant_marker_pvals)), 2)
malignant_marker_lfcs <- lapply(de_malignant_timepoint_files, function(f) {
  de_res <- readRDS(f)
  pvals <- (de_res$malignant$gene %>% 
              dplyr::filter(Symbol %in% malignant_markers))$logFC
  return(pvals)
})
min_malignant_lfc <- signif(min(unlist(malignant_marker_lfcs)), 2)

## Cell proportions

message("FL cell proportion stats ...")

b_rel_proportions <- plyr::rbind.fill(lapply(unique(sce_follicular_b$patient), function(pat) {
  sce_subset <- sce_follicular_b %>%
    scater::filter(patient == pat)
  tab <- with(colData(sce_subset), table(timepoint, malignant_status_manual))
  return(data.frame(patient=pat,
                    as.data.frame(tab/rowSums(tab)) %>%
                      dplyr::rename(prop=Freq) %>%
                      dplyr::mutate(percentage=round(prop*100, 1))))
}))

T1_nonmalignant_b_res <- b_rel_proportions %>%
  dplyr::filter(timepoint == "T1",
                malignant_status_manual == "nonmalignant")
T1_nonmalignant_b_pct <- as.list(T1_nonmalignant_b_res$percentage)
names(T1_nonmalignant_b_pct) <- T1_nonmalignant_b_res$patient

T2_nonmalignant_b_res <- b_rel_proportions %>%
  dplyr::filter(timepoint == "T2",
                malignant_status_manual == "nonmalignant")
T2_nonmalignant_b_pct <- as.list(T2_nonmalignant_b_res$percentage)
names(T2_nonmalignant_b_pct) <- T2_nonmalignant_b_res$patient

## DE significance

message("FL DE stats ...")

interesting_genes <- c(
  "CD69",
  "IFNG",
  "GZMA",
  "PRF1"
)

interesting_significance <- plyr::rbind.fill(lapply(de_timepoint_files, function(f) {
  res <- readRDS(f)
  celltype <- basename(dirname(f))
  patient <- basename(tools::file_path_sans_ext(f))
  df <- res$gene %>% 
    dplyr::filter(Symbol %in% interesting_genes)
  return(data.frame(patient=patient, celltype=celltype, df))
}))

df_to_list <- function(df) {
  x <- as.list(df$adj.P.Val)
  names(x) <- df$Symbol
  return(x)
}

fl1018_cytotoxic_pvals <- df_to_list(interesting_significance %>% 
                                       subset(patient == "FL1018" &
                                                celltype == "cytotoxic"))
fl1018_tfh_pvals <- df_to_list(interesting_significance %>% 
                                 subset(patient == "FL1018" &
                                          celltype == "follicular_helper"))
fl1018_cd4_pvals <- df_to_list(interesting_significance %>% 
                                 subset(patient == "FL1018" &
                                          celltype == "helper"))
fl1018_malignant_pvals <- df_to_list(interesting_significance %>% 
                                       subset(patient == "FL1018" &
                                                celltype == "malignant"))
fl2001_malignant_pvals <- df_to_list(interesting_significance %>% 
                                       subset(patient == "FL2001" &
                                                celltype == "malignant"))

max_t_cd69_pval <- signif(max(fl1018_tfh_pvals$CD69, fl1018_cd4_pvals$CD69), 2)



fgsea_res <- lapply(de_timepoint_fgsea_files, function(f) {
  readRDS(f)
})
names(fgsea_res) <- tools::file_path_sans_ext(basename(de_timepoint_fgsea_files))

# FL2001_E2F_padj <- (fgsea_res$FL2001$pathway %>% 
#   dplyr::filter(str_detect(pathway, "E2F")))$padj %>%
#   signif(2)

fl1018_prolif_replicative_pvals <- max((fgsea_res$FL1018$pathway %>%
                                          dplyr::filter(pathway %in% 
                                                          c("HALLMARK_E2F_TARGETS",
                                                            "HALLMARK_G2M_CHECKPOINT")))$padj %>%
                                         signif(2))

fl1018_prolif_replicative_nes <- min((fgsea_res$FL1018$pathway %>%
                                          dplyr::filter(pathway %in% 
                                                          c("HALLMARK_E2F_TARGETS",
                                                            "HALLMARK_G2M_CHECKPOINT")))$NES %>%
                                         signif(2))

fl2001_mitotic_spindle_padj <- (fgsea_res$FL2001$pathway %>%
  dplyr::filter(str_detect(pathway, "MITOTIC_SPINDLE")))$padj %>%
  signif(2)

fl2001_mitotic_spindle_nes <- (fgsea_res$FL2001$pathway %>%
                                  dplyr::filter(str_detect(pathway, "MITOTIC_SPINDLE")))$NES %>%
  signif(3)

fl2001_myc_vone_padj <- (fgsea_res$FL2001$pathway %>%
                           dplyr::filter(str_detect(pathway, "MYC_TARGETS_V1")))$padj %>%
  signif(2)

fl2001_etwof_gtwom_padj <- min((fgsea_res$FL2001$pathway %>%
                                  dplyr::filter(pathway %in% 
                                                  c("HALLMARK_E2F_TARGETS",
                                                    "HALLMARK_G2M_CHECKPOINT")))$padj %>%
                                 signif(2))

## Replication stats

fl1018_prolif_pvals <- fgsea_res$FL1018$gene %>% 
  dplyr::filter(Symbol %in% c("MKI67", "TOP2A"), 
                contrast == "T2") %>%
  dplyr::rename(adj.P.Val=FDR) %>%
  dplyr::mutate(adj.P.Val=signif(adj.P.Val,2)) %>%
  df_to_list()

fl2001_prolif_pvals <- fgsea_res$FL2001$gene %>% 
  dplyr::filter(Symbol %in% c("MKI67", "TOP2A"), 
                contrast == "T2") %>%
  dplyr::rename(adj.P.Val=FDR) %>%
  dplyr::mutate(adj.P.Val=signif(adj.P.Val,2)) %>%
  df_to_list()



## Timepoint variability

timepoint_correlations <- plyr::rbind.fill(lapply(unique(sce_follicular$patient), function(pat) {
  sce_follicular_subset <- sce_follicular %>% 
    scater::filter(patient == pat)
  
  corsummary <- plyr::rbind.fill(lapply(unique(sce_follicular_subset$malignant_status_manual), function(status) {
    sce_tmp1 <- sce_follicular_subset %>%
      scater::filter(malignant_status_manual == status,
                     timepoint == "T1")
    sce_tmp2 <- sce_follicular_subset %>%
      scater::filter(malignant_status_manual == status,
                     timepoint == "T2")
    
    summary <- data.frame(cor=cor(rowMeans(logcounts(sce_tmp1)), rowMeans(logcounts(sce_tmp2))),
               malignant_status_manual=status)
    
    return(summary)
  }))
  return(data.frame(patient=pat, corsummary))
}))

fl1018_malignant_timepoint_var <- (timepoint_correlations %>%
                                     dplyr::filter(malignant_status_manual == "malignant",
                                                   patient == "FL1018"))$cor %>%
  signif(3)
fl2001_malignant_timepoint_var <- (timepoint_correlations %>%
                                     dplyr::filter(malignant_status_manual == "malignant",
                                                   patient == "FL2001"))$cor %>% 
  signif(3)

hla_genes <- c("HLA-A", "HLA-B",
               "HLA-C", "B2M",
               "HLA-DRA", "HLA-DRB1")

de_malignant_timepoint_res <- lapply(de_malignant_timepoint_files, function(f) {
  readRDS(f)
})
names(de_malignant_timepoint_res) <- tools::file_path_sans_ext(basename(de_malignant_timepoint_files))

fl1018_hla_malignant_pvals <- de_malignant_timepoint_res$FL1018$malignant$gene %>%
  dplyr::filter(Symbol %in% hla_genes)

fl2001_hla_malignant_pvals <- de_malignant_timepoint_res$FL2001$malignant$gene %>%
  dplyr::filter(Symbol %in% hla_genes)

max_hla_malignant_pval <- max(fl1018_hla_malignant_pvals$FDR, fl2001_hla_malignant_pvals$FDR) %>%
  signif(2)

fl1018_malignant <- readRDS(de_timepoint_files[basename(dirname(de_timepoint_files)) == "malignant" & str_detect(de_timepoint_files, "FL1018")])

fl1018_class1_MHC_p <- (fl1018_malignant$pathway$down@result %>%
  dplyr::filter(str_detect(Description, "Antigen Presentation: Folding, assembly and peptide loading of class I MHC")))$p.adjust %>%
  signif(2)

fl1018_hla_timepoint_maxpval <- (fl1018_malignant$gene %>%
                       dplyr::filter(Symbol %in% hla_genes))$FDR %>%
  max %>%
  signif(2)

fl1018_hla_timepoint_minlfc <- (fl1018_malignant$gene %>%
                                   dplyr::filter(Symbol %in% hla_genes))$logFC %>%
  max %>%
  signif(2)

## HGSC

message("HGSC stats ...")

remap_hgsc_names <- function(x) {
  x %>% plyr::mapvalues(c("Left Ovary", "Right Ovary"),
                        c("leftovary", "rightovary"))
}

hgsc_total_cell_count <- ncol(sce_hgsc_raw)
hgsc_cell_counts <- as.list(table(sce_hgsc_raw$patient))
hgsc_cell_counts_by_sample <- as.list(table(sce_hgsc_raw$dataset %>%
                                              remap_hgsc_names()))

hgsc_num_celltypes <- length(setdiff(unique(sce_hgsc$celltype), "other"))

hgsc_celltype_table <- with(colData(sce_hgsc), table(dataset, celltype))
hgsc_celltype_props <- hgsc_celltype_table/rowSums(hgsc_celltype_table)
immune_prevalence <- (rowSums(hgsc_celltype_props[,c("B cells", "T cells", "Monocyte/Macrophage")]) * 100) %>%
  signif(digits = 2)
names(immune_prevalence) <- names(immune_prevalence) %>% 
  remap_hgsc_names()
myeloid_prevalence <- (hgsc_celltype_props[,c("Monocyte/Macrophage")]* 100) %>%
  signif(digits = 2)
names(myeloid_prevalence) <- names(myeloid_prevalence) %>% 
  remap_hgsc_names()

myeloid_immune_prop <- round(hgsc_celltype_props[,c("Monocyte/Macrophage")]/rowSums(hgsc_celltype_props[,c("B cells", "T cells", "Monocyte/Macrophage")])*100, 1)
names(myeloid_immune_prop) <- names(myeloid_immune_prop) %>% 
  remap_hgsc_names()

de_site_epithelial_fgsea_res <- readRDS(Sys.glob(file.path(de_site_fgsea_dir, "epithelial", "*.rds")))
epithelial_emt_padj <- (de_site_epithelial_fgsea_res$pathway %>%
  dplyr::filter(str_detect(pathway, "MESENCHYMAL")))$padj %>%
  signif(2)

epithelial_emt_nes <- ((de_site_epithelial_fgsea_res$pathway %>%
                          dplyr::filter(str_detect(pathway, "MESENCHYMAL")))$NES * -1) %>%
  signif(3)

epithelial_ecad_padj <- (de_site_epithelial_fgsea_res$gene %>%
                           dplyr::filter(Symbol == "CDH1"))$FDR[1] %>%
  signif(2)

epithelial_ecad_lfc <- (de_site_epithelial_fgsea_res$gene %>%
                           dplyr::filter(Symbol == "CDH1"))$logFC.Right.ovary[1] %>%
  signif(2)

de_cluster_epithelial_fgsea_res <- readRDS(Sys.glob(file.path(de_epithelial_dir, "*.rds")))

hypoxia_padj <- (de_cluster_epithelial_fgsea_res$pathway %>%
  dplyr::filter(str_detect(pathway, "HYPOXIA"),
                clust1 == "2",
                clust2 %in% c("0", "4")))$padj %>%
  signif(2) %>%
  max

hypoxia_nes <- (de_cluster_epithelial_fgsea_res$pathway %>%
                   dplyr::filter(str_detect(pathway, "HYPOXIA"),
                                 clust1 == "2",
                                 clust2 %in% c("0", "4")))$NES %>%
  signif(3) %>%
  min

apoptosis_padj <- (de_cluster_epithelial_fgsea_res$pathway %>%
                     dplyr::filter(str_detect(pathway, "APOPTOSIS"),
                                   clust1 == "2",
                                   clust2 %in% c("0", "4")))$padj %>%
  signif(2) %>%
  max

glycolysis_padj <- (de_cluster_epithelial_fgsea_res$pathway %>%
                     dplyr::filter(str_detect(pathway, "GLYCOLYSIS"),
                                   clust1 == "2",
                                   clust2 %in% c("0", "4")))$padj %>%
  signif(2) %>%
  max

## Koh stats

sce_koh <- readRDS(koh_annotated_path)
categorical_palettes <- cat_palettes()
factor_orderings <- factor_orders()

evaluation_measures <- sce_koh@metadata$evaluation_measures %>%
  dplyr::mutate(full_label=ifelse(is.na(gene_set),
                                  as.character(clustering_method),
                                  paste0(clustering_method, "_", gene_set)))

koh_external_methods <- c("SC3_full",
                          "seurat_0.8_full",
                          "seurat_0.8_markers",
                          "seurat_1.2_full",
                          "seurat_1.2_markers")

external_stats <- evaluation_measures %>%
  dplyr::filter(full_label %in% koh_external_methods)

koh_scina_conttable <- table(sce_koh$celltype, scina_fit_koh$cell_labels)
koh_scina_accuracy <- (sum(diag(koh_scina_conttable))/sum(koh_scina_conttable)) %>%
  round(3)
koh_scina_macrof1 <- macroF1(koh_scina_conttable) %>%
  round(3)

koh_external_best_macrof1 <- c(external_stats$macro_f1, koh_scina_macrof1) %>% max %>%
  round(3)
koh_external_best_accuracy <- c(external_stats$accuracy, koh_scina_accuracy) %>% max %>%
  round(3)

koh_cellassign_conttable <- table(sce_koh$celltype, cellassign_fit_koh$cell_type)
koh_cellassign_accuracy <- (sum(diag(koh_cellassign_conttable))/sum(koh_cellassign_conttable)) %>%
  round(3)
koh_cellassign_macrof1 <- macroF1(koh_cellassign_conttable) %>%
  round(3)

# koh_cellassign_macrof1 <- evaluation_measures[evaluation_measures$clustering_method == "cellassign-shrinkage",]$macro_f1 %>%
#   round(3)
# koh_cellassign_accuracy <- evaluation_measures[evaluation_measures$clustering_method == "cellassign-shrinkage",]$accuracy %>%
#   round(3)

sce_koh_aps_mps <- sce_koh %>%
  scater::filter(celltype %in% c("APS", "MPS"))
aps_mps_correct <- length(which(sce_koh_aps_mps$celltype == sce_koh_aps_mps$cellassign_cluster))
aps_mps_total <- ncol(sce_koh_aps_mps)


## Follicular cytotoxic cell numbers

cytotoxic_cell_counts <- table(sce_follicular_cytotoxic$dataset)

## Follicular malignant/nonmalignant B numbers

sce_follicular_b_filtered <- sce_follicular_b %>%
  scater::filter(celltype_full %in% c("B cells", "B cells (malignant)"))

malignant_counts <- table((sce_follicular_b_filtered %>% scater::filter(celltype_full == "B cells (malignant)"))$dataset)
nonmalignant_counts <- table((sce_follicular_b_filtered %>% scater::filter(celltype_full == "B cells"))$dataset)

## Liver

message("Liver stats ...")

liver_cell_count <- ncol(sce_liver)
liver_patient_count <- length(unique(sce_liver$patient))
liver_celltype_count <- length(unique(sce_liver$celltype))

liver_fit_stats <- lapply(names(liver_fits), function(x) {
  fit <- liver_fits[[x]]
  if (class(fit) == "cellassign_fit") {
    celltypes <- setdiff(colnames(fit$mle_params$gamma), "other")
    assignments <- fit$cell_type %>% plyr::mapvalues(from = c("other"), to = c("unknown"))
  } else {
    celltypes <- setdiff(rownames(fit$probabilities), "unknown")
    assignments <- fit$cell_labels
  }
  
  sce_subset <- sce_liver %>%
    scater::filter(celltype %in% liver_data_celltypes[[x]])
  
  sce_subset$celltype[!sce_subset$celltype %in% celltypes] <- "unknown"
  
  celltype_levels <- union(unique(sce_subset$celltype), assignments)
  
  cont_table <- table(factor(sce_subset$celltype, levels = celltype_levels), 
                      factor(assignments, levels = celltype_levels))
  acc <- signif(sum(diag(cont_table))/sum(cont_table) * 100, 3)
  micro_f1 <- signif(microF1(cont_table), 3)
  
  return(list('correct'=as.integer(sum(diag(cont_table))), 'total'=as.integer(sum(cont_table)), 'acc'=acc, 'fone'=micro_f1))
})
names(liver_fit_stats) <- names(liver_fits)

liver_fit_stats_flat <- unlist(liver_fit_stats)
names(liver_fit_stats_flat) <- str_replace_all(names(liver_fit_stats_flat), "[_\\.]", "")
liver_fit_stats_final <- as.list(liver_fit_stats_flat)

## Combined

message("HGSC + FL stats ...")

hgsc_fl_combination_celltypes <- data.frame(
  Sample=sce_fl_hgsc_merged$Sample,
  Barcode=sce_fl_hgsc_merged$Barcode,
  combined=combination_fits$cellassign_fl_hgsc$cell_type,
  stringsAsFactors = FALSE
)

fl_celltypes <- data.frame(
  Sample=sce_follicular$Sample,
  Barcode=sce_follicular$Barcode,
  separate=combination_fits$cellassign_fl_combination$cell_type,
  type='follicular',
  stringsAsFactors = FALSE
)

hgsc_celltypes <- data.frame(
  Sample=sce_hgsc$Sample,
  Barcode=sce_hgsc$Barcode,
  separate=combination_fits$cellassign_hgsc_combination$cell_type,
  type='hgsc',
  stringsAsFactors = FALSE
)

hgsc_fl_combination_merged <- hgsc_fl_combination_celltypes %>%
  dplyr::inner_join(rbind(fl_celltypes, hgsc_celltypes))

hgsc_fl_combination_merged_summary <- hgsc_fl_combination_merged %>%
  dplyr::mutate(combined=ifelse(type == "hgsc", plyr::mapvalues(combined,
                                                                c("B cells (kappa)", "B cells (lambda)"),
                                                                rep("B cells", 2)), combined)) %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(correct=sum(combined == separate),
                   total=length(combined),
                   accuracy=signif(sum(combined == separate)/length(combined) * 100, 3),
                   micro_f1=signif(microF1(table(factor(combined, levels = unique(c(combined, separate))), 
                                                 factor(separate, levels = unique(c(combined, separate))))), 3))

combination_stats <- list(
  'hgscflcells'=nrow(hgsc_fl_combination_celltypes),
  'hgscflflcorrect'=subset(hgsc_fl_combination_merged_summary, type == "follicular")$correct,
  'hgscflfltotal'=subset(hgsc_fl_combination_merged_summary, type == "follicular")$total,
  'hgscflflacc'=subset(hgsc_fl_combination_merged_summary, type == "follicular")$accuracy,
  'hgscflflfone'=subset(hgsc_fl_combination_merged_summary, type == "follicular")$micro_f1,
  'hgscflhgsccorrect'=subset(hgsc_fl_combination_merged_summary, type == "hgsc")$correct,
  'hgscflhgsctotal'=subset(hgsc_fl_combination_merged_summary, type == "hgsc")$total,
  'hgscflhgscacc'=subset(hgsc_fl_combination_merged_summary, type == "hgsc")$accuracy,
  'hgscflhgscfone'=subset(hgsc_fl_combination_merged_summary, type == "hgsc")$micro_f1
)

## Hierarchical

message("Hierarchical stats ...")

fl_noother_celltypes <- data.frame(
  Sample=sce_follicular$Sample,
  Barcode=sce_follicular$Barcode,
  combined=combination_fits$cellassign_fl_noother$cell_type,
  stringsAsFactors = FALSE
)

fl_top_celltypes <- data.frame(
  Sample=sce_follicular$Sample,
  Barcode=sce_follicular$Barcode,
  hierarchical=combination_fits$cellassign_fl_top$cell_type,
  stringsAsFactors = FALSE
)

fl_bottom_celltypes <- data.frame(
  Sample=sce_follicular_nonmalignant_b$Sample,
  Barcode=sce_follicular_nonmalignant_b$Barcode,
  hierarchical=plyr::mapvalues(combination_fits$cellassign_fl_bottom$cell_type,
                               from = c("IGKC", "IGLC"),
                               to = c("B cells (kappa)", "B cells (lambda)")),
  stringsAsFactors = FALSE
)

fl_bottom_combined <- fl_noother_celltypes %>%
  dplyr::inner_join(fl_bottom_celltypes)

fl_top_combined <- fl_noother_celltypes %>%
  dplyr::mutate(combined=plyr::mapvalues(combined, from = c("B cells (kappa)", "B cells (lambda)"),
                                         to = rep("B cells", 2))) %>%
  dplyr::inner_join(fl_top_celltypes)

fl_bottom_conttable <- with(fl_bottom_combined, table(combined, hierarchical))
fl_bottom_correct <- sum(diag(fl_bottom_conttable))
fl_bottom_total <- sum(fl_bottom_conttable)
fl_bottom_accuracy <- signif(fl_bottom_correct/fl_bottom_total * 100, 3)
fl_bottom_f1 <- signif(microF1(fl_bottom_conttable), 3)

fl_top_conttable <- with(fl_top_combined, table(combined, hierarchical))
fl_top_correct <- sum(diag(fl_top_conttable))
fl_top_total <- sum(fl_top_conttable)
fl_top_accuracy <- signif(fl_top_correct/fl_top_total * 100, 3)
fl_top_f1 <- signif(microF1(fl_top_conttable), 3)

sce_koh$Sample <- "Temp"
sce_koh$Barcode <- paste0("Barcode", 1:ncol(sce_koh))

koh_noother_celltypes <- data.frame(
  Sample=sce_koh$Sample,
  Barcode=sce_koh$Barcode,
  combined=combination_fits$cellassign_koh_noother$cell_type,
  stringsAsFactors = FALSE
)

koh_top_celltypes <- data.frame(
  Sample=sce_koh$Sample,
  Barcode=sce_koh$Barcode,
  hierarchical=combination_fits$cellassign_koh_top$cell_type,
  stringsAsFactors = FALSE
)

sce_koh_aps_mps <- sce_koh[,combination_fits$cellassign_koh_top$cell_type == "APS_MPS"]

koh_bottom_celltypes <- data.frame(
  Sample=sce_koh_aps_mps$Sample,
  Barcode=sce_koh_aps_mps$Barcode,
  hierarchical=combination_fits$cellassign_koh_bottom$cell_type,
  stringsAsFactors = FALSE
)

koh_bottom_combined <- koh_noother_celltypes %>%
  dplyr::inner_join(koh_bottom_celltypes)

koh_top_combined <- koh_noother_celltypes %>%
  dplyr::mutate(combined=plyr::mapvalues(combined, from = c("APS", "MPS"),
                                         to = rep("APS_MPS", 2))) %>%
  dplyr::inner_join(koh_top_celltypes)

koh_bottom_conttable <- with(koh_bottom_combined, table(combined, hierarchical))
koh_bottom_correct <- sum(diag(koh_bottom_conttable))
koh_bottom_total <- sum(koh_bottom_conttable)
koh_bottom_accuracy <- signif(koh_bottom_correct/koh_bottom_total * 100, 3)
koh_bottom_f1 <- signif(microF1(koh_bottom_conttable), 3)

koh_top_conttable <- with(koh_top_combined, table(combined, hierarchical))
koh_top_correct <- sum(diag(koh_top_conttable))
koh_top_total <- sum(koh_top_conttable)
koh_top_accuracy <- signif(koh_top_correct/koh_top_total * 100, 3)
koh_top_f1 <- signif(microF1(koh_top_conttable), 3)

hierarchical_stats <- list(
  'flbottomcorrect'=fl_bottom_correct,
  'fltopcorrect'=fl_top_correct,
  'flbottomtotal'=fl_bottom_total,
  'fltoptotal'=fl_top_total,
  'flbottomfone'=fl_bottom_f1,
  'fltopfone'=fl_top_f1,
  'flbottomacc'=fl_bottom_accuracy,
  'fltopacc'=fl_top_accuracy,
  'kohbottomcorrect'=koh_bottom_correct,
  'kohtopcorrect'=koh_top_correct,
  'kohbottomtotal'=koh_bottom_total,
  'kohtoptotal'=koh_top_total,
  'kohbottomfone'=koh_bottom_f1,
  'kohtopfone'=koh_top_f1,
  'kohbottomacc'=koh_bottom_accuracy,
  'kohtopacc'=koh_top_accuracy
)

## CellBench

message("CellBench stats ...")

pure_fits <- cellbench_fits[str_detect(names(cellbench_fits), "pure")]

pure_fit_stats <- lapply(names(pure_fits), function(x) {
  fit <- pure_fits[[x]]
  sce <- cellbench_sets[[x]]
  
  cont_table <- table(sce$cell_line,
                      fit$cell_type)
  acc <- signif(sum(diag(cont_table))/sum(cont_table) * 100, 3)
  micro_f1 <- signif(microF1(cont_table), 3)
  
  return(list('correct'=as.integer(sum(diag(cont_table))), 'total'=as.integer(sum(cont_table)), 'acc'=acc, 'fone'=micro_f1))
})
names(pure_fit_stats) <- names(pure_fits)

pure_fit_stats_flat <- unlist(pure_fit_stats)
names(pure_fit_stats_flat) <- str_replace_all(names(pure_fit_stats_flat), "[_\\.]", "")
pure_fit_stats_final <- as.list(pure_fit_stats_flat)

## Stats

stats <- list(
  minDeltaR=min_delta_rval,
  kohMarkers=num_koh_markers,
  flTotalCellCount=fl_total_cell_count,
  flCellCountsByPatient=fl_cell_counts,
  flCellCountsBySample=fl_cell_counts_by_sample,
  lkNonOther=lk_nonother,
  lkKappa=lk_kappa,
  lkLambda=lk_lambda,
  lkKappaPct=lk_kappa_pct,
  lkKappaPatient=lk_kappa_patient,
  lkLambdaPatient=lk_lambda_patient,
  lkNonOtherPatient=lk_nonother_patient,
  lkKappaPctPatient=lk_kappa_pct_patient,
  maxMalignantPval=max_malignant_pval,
  negminmalignantlfc=-1*min_malignant_lfc,
  fl1018MalignantVar=fl1018_malignant_timepoint_var,
  fl2001MalignantVar=fl2001_malignant_timepoint_var,
  T1NonMalignantB=T1_nonmalignant_b_pct,
  T2NonMalignantB=T2_nonmalignant_b_pct,
  maxCDsixninePval=max_t_cd69_pval,
  fl1018Prolifp=fl1018_prolif_pvals,
  fl2001Prolifp=fl2001_prolif_pvals,
  maxHLAMalignantPval=max_hla_malignant_pval,
  fl1018classoneMHCp=fl1018_class1_MHC_p,
  fl1018maxHLAtimepointPval=fl1018_hla_timepoint_maxpval,
  fl1018minhlatimepointlfc=fl1018_hla_timepoint_minlfc,
  fl2001mitoticspindlePval=fl2001_mitotic_spindle_padj,
  fl2001mitoticspindlenes=fl2001_mitotic_spindle_nes,
  fl1018prolifreplicativepvals=fl1018_prolif_replicative_pvals,
  fl1018prolifreplicativenes=fl1018_prolif_replicative_nes,
  fl2001myconepval=fl2001_myc_vone_padj,
  fl2001etwofgtwompadj=fl2001_etwof_gtwom_padj,
  hgscTotalCellCount=hgsc_total_cell_count,
  hgscCellCountsByPatient=hgsc_cell_counts,
  hgscCellCountsBySample=hgsc_cell_counts_by_sample,
  hgscNumCelltypes=hgsc_num_celltypes,
  hgscImmunePrevalence=as.list(immune_prevalence),
  myeloidImmuneProp=as.list(myeloid_immune_prop),
  epithelialemtpval=epithelial_emt_padj,
  epithelialemtnes=epithelial_emt_nes,
  epithelialecadpval=epithelial_ecad_padj,
  epithelialecadlfc=epithelial_ecad_lfc,
  epithelialclusterhypoxiapval=hypoxia_padj,
  epithelialclusterhypoxianes=hypoxia_nes,
  epithelialclusterapoptosisglycolysispval=max(glycolysis_padj, apoptosis_padj),
  kohexternalbestfone=koh_external_best_macrof1,
  kohexternalbestacc=koh_external_best_accuracy,
  kohcellassignfone=koh_cellassign_macrof1,
  kohcellassignacc=koh_cellassign_accuracy,
  kohapsmpstotal=aps_mps_total,
  kohapsmpscorrect=aps_mps_correct,
  flcytotoxiccellcounts=cytotoxic_cell_counts,
  flmalignanttimepointcounts=malignant_counts,
  flnonmalignanttimepointcounts=nonmalignant_counts,
  nlivercells=liver_cell_count,
  nlivercelltypes=liver_celltype_count,
  nliverpatients=liver_patient_count,
  nliverpanglaodbmarkers=length(unique(unlist(liver_panglaodb_markers)))
)

## Combine sub-stat data frames
stats <- c(stats, liver_fit_stats_final, combination_stats, hierarchical_stats, pure_fit_stats_final)

# Write outputs

stats_flat <- as.list(unlist(stats))
names(stats_flat) <- names(stats_flat) %>%
  tolower() %>%
  stringr::str_replace("\\.", "") %>%
  stringr::str_replace("fl1018", "fltrans") %>%
  stringr::str_replace("fl2001", "flprog") %>%
  stringr::str_replace("t2", "two") %>%
  stringr::str_replace("t1", "one") %>%
  stringr::str_replace("mki67", "mki") %>%
  stringr::str_replace("top2a", "topo") %>%
  stringr::str_replace("voa11543", "hgscpatientone")

list_to_tex <- function(x) {
  strs <- lapply(names(x), function(tag) {
    ## Remove underscores because this will screw up TeX
    tag_p <- stringr::str_replace_all(tag, "_", "")
    command <- paste("\\newcommand{\\", tag_p, "}{", x[[tag]], "}", sep="")
    return(command)
  })
  ## Store in dataframe because of newline ambiguity
  final_str <- data.frame(texstr=unlist(strs))
  return(final_str)
}

output_string <- list_to_tex(stats_flat)

write.table(output_string, file=args$outfname, row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Completed.\n")

