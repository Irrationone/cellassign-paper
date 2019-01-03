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
parser$add_argument('--koh_annotated', metavar='FILE', type='character',
                    help="Annotated SCE of Koh")
parser$add_argument('--koh_rho', type='character', metavar='FILE',
                    help="Koh et al. rho matrix")
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
koh_annotated_path <- args$koh_annotated

# HGSC data

sce_hgsc <- readRDS(args$sce_hgsc)
sce_hgsc_raw <- readRDS(args$sce_hgsc_raw)
de_site_dir <- args$de_site
de_site_fgsea_dir <- args$de_site_fgsea
de_epithelial_dir <- args$de_epithelial

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

## HGSC

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

epithelial_ecad_padj <- (de_site_epithelial_fgsea_res$gene %>%
                           dplyr::filter(Symbol == "CDH1"))$FDR[1] %>%
  signif(2)

de_cluster_epithelial_fgsea_res <- readRDS(Sys.glob(file.path(de_epithelial_dir, "*.rds")))

hypoxia_padj <- (de_cluster_epithelial_fgsea_res$pathway %>%
  dplyr::filter(str_detect(pathway, "HYPOXIA"),
                clust1 == "2",
                clust2 %in% c("0", "4")))$padj %>%
  signif(2) %>%
  max

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

koh_external_best_macrof1 <- external_stats$macro_f1 %>% max %>%
  round(3)
koh_external_best_accuracy <- external_stats$accuracy %>% max %>%
  round(3)

koh_cellassign_macrof1 <- evaluation_measures[evaluation_measures$clustering_method == "cellassign-shrinkage",]$macro_f1 %>%
  round(3)
koh_cellassign_accuracy <- evaluation_measures[evaluation_measures$clustering_method == "cellassign-shrinkage",]$accuracy %>%
  round(3)

sce_koh_aps_mps <- sce_koh %>%
  scater::filter(celltype %in% c("APS", "MPS"))
aps_mps_correct <- length(which(sce_koh_aps_mps$celltype == sce_koh_aps_mps$cellassign_cluster))
aps_mps_total <- ncol(sce_koh_aps_mps)


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
  fl2001etwofgtwompadj=fl2001_etwof_gtwom_padj,
  hgscTotalCellCount=hgsc_total_cell_count,
  hgscCellCountsByPatient=hgsc_cell_counts,
  hgscCellCountsBySample=hgsc_cell_counts_by_sample,
  hgscNumCelltypes=hgsc_num_celltypes,
  hgscImmunePrevalence=as.list(immune_prevalence),
  myeloidImmuneProp=as.list(myeloid_immune_prop),
  epithelialemtpval=epithelial_emt_padj,
  epithelialecadpval=epithelial_ecad_padj,
  epithelialclusterhypoxiapval=hypoxia_padj,
  epithelialclusterapoptosisglycolysispval=max(glycolysis_padj, apoptosis_padj),
  kohexternalbestfone=koh_external_best_macrof1,
  kohexternalbestacc=koh_external_best_accuracy,
  kohcellassignfone=koh_cellassign_macrof1,
  kohcellassignacc=koh_cellassign_accuracy,
  kohapsmpstotal=aps_mps_total,
  kohapsmpscorrect=aps_mps_correct
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

