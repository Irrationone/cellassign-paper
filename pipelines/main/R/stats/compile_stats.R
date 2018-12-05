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

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create compiled stats file.")
parser$add_argument('--sce_follicular', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for FL")
parser$add_argument('--sce_follicular_bcells', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for FL")
parser$add_argument('--sce_follicular_raw', metavar='FILE', type='character',
                    help="Path to raw SingleCellExperiment RDS")
parser$add_argument('--patient_progression', metavar='FILE', type='character',
                    help="Patient events")
parser$add_argument('--deprob_result_dir', metavar='DIR', type='character',
                    help="Path to simulation result directory for DE prob")
parser$add_argument('--delta_deprobs', type='double', nargs ='+',
                    help="DE probs to show delta plots for.")
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
parser$add_argument('--koh_rho', type='character', metavar='FILE',
                    help="Koh et al. rho matrix")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for stats file")
args <- parser$parse_args()

sce_follicular <- readRDS(args$sce_follicular)
sce_follicular_raw <- readRDS(args$sce_follicular_raw)
sce_follicular_nonmalignant_b <- readRDS(args$sce_follicular_nonmalignant_bcell)
sce_follicular_b <- readRDS(args$sce_follicular_bcells)
deprob_result_dir <- args$deprob_result_dir
delta_deprobs <- unlist(args$delta_deprobs)
cellassign_lambda_kappa_results <- readRDS(args$cellassign_lambda_kappa)

de_timepoint_dir <- args$de_timepoint_dir
de_timepoint_fgsea_dir <- args$de_timepoint_fgsea_dir
de_malignant_timepoint_dir <- args$de_malignant_timepoint_dir

de_timepoint_files <- Sys.glob(file.path(de_timepoint_dir, "*", "*"))
de_timepoint_fgsea_files <- Sys.glob(file.path(de_timepoint_fgsea_dir, "malignant", "*"))
de_malignant_timepoint_files <- Sys.glob(file.path(de_malignant_timepoint_dir, "*", "*"))

koh_rho <- read.table(args$koh_rho, sep = "\t", header = TRUE, row.names = 1) %>%
  tibble::rownames_to_column(var = "Gene")

## Delta stats

de_eval_measures <- load_annotation_files(deprob_result_dir, pattern = "*_eval_measures.tsv")
de_deltas <- load_annotation_files(deprob_result_dir, pattern = "*_delta_compare.tsv")

delta_table <- de_deltas %>% 
  dplyr::filter(de_prob %in% delta_deprobs) 

rvals <- compute_pvals_subsets(delta_table,
                               facet_vars = c("de_prob", "clustering_method"),
                               formula = ~ true_delta + inferred_delta,
                               corfun = cor.test,
                               output = "estimate")

min_delta_rval <- round(min(rvals$estimate), 3)

## Koh rho

num_koh_markers <- nrow(koh_rho)

## FL cell counts

fl_total_cell_count <- ncol(sce_follicular_raw)
fl_cell_counts <- as.list(table(sce_follicular_raw$patient))
fl_cell_counts_by_sample <- as.list(table(sce_follicular_raw$dataset))

## Lambda/kappa stats

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

malignant_markers <- c("BCL2", "BCL6")
malignant_marker_pvals <- lapply(de_malignant_timepoint_files, function(f) {
  de_res <- readRDS(f)
  pvals <- (de_res$malignant$gene %>% 
              dplyr::filter(Symbol %in% malignant_markers))$adj.P.Val
  return(pvals)
})
max_malignant_pval <- signif(max(unlist(malignant_marker_pvals)), 2)

## Cell proportions

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
                                                          c("HALLMARK_MYC_TARGETS_V1",
                                                            "HALLMARK_MYC_TARGETS_V2",
                                                            "HALLMARK_E2F_TARGETS",
                                                            "HALLMARK_G2M_CHECKPOINT")))$padj %>%
                                         signif(2))

fl2001_mitotic_spindle_padj <- (fgsea_res$FL2001$pathway %>%
  dplyr::filter(str_detect(pathway, "MITOTIC_SPINDLE")))$padj %>%
  signif(2)

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
  fl2001mitoticspindlePval=fl2001_mitotic_spindle_padj,
  fl1018prolifreplicativepvals=fl1018_prolif_replicative_pvals,
  fl2001myconepval=fl2001_myc_vone_padj,
  fl2001etwofgtwompadj=fl2001_etwof_gtwom_padj
)

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
  stringr::str_replace("top2a", "topo")

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

