# Figure for proportion validation

# Duo validation
# HGSC validation

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(cowplot)
library(pheatmap)
library(Matrix)
library(DuoClustering2018)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Create proportion validation figure.")
parser$add_argument('--koh_bulkrna', type = 'character', metavar = 'FILE',
                    help="Bulk RNA file for Koh")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for PDF plot")
args <- parser$parse_args()

# Koh analysis

koh_bulk_scrna_celltypes <- c("hESC", "APS", "MPS", "DLL1pPXM", "ESMT", "Sclrtm", "D5CntrlDrmmtm",
                              "D2LtM")

koh_common_celltype_map <- c(
  'D0 H7 hESC'='hESC', 
  'D1 APS'='APS', 
  'D1 MPS'='MPS', 
  'D2 DLL1pos PXM'="DLL1pPXM", 
  'D3 Somite'="ESMT", 
  'D6 Sclerotome'="Sclrtm", 
  'D5 Dermomyotome'="D5CntrlDrmmtm",
  'D2 LatM'="D2LtM"
)

sce_koh <- DuoClustering2018::sce_filteredExpr10_Koh()

sce_koh <- sce_koh %>%
  scater::mutate(celltype = str_replace_all(phenoid, "^H7(_derived_|_dreived_)?", "")) 

sce_koh_filtered <- sce_koh %>% 
  scater::filter(celltype %in% koh_bulk_scrna_celltypes)

sce_koh_normalized <- normalize(sce_koh_filtered)
sce_koh_normalized <- runPCA(sce_koh_normalized, ntop = 500, 
                             ncomponents = 50, exprs_values = "logcounts")
sce_koh_normalized <- runTSNE(sce_koh_normalized, use_dimred = "PCA", 
                              n_dimred = 50, ncomponents = 2)
sce_koh_normalized <- runUMAP(sce_koh_normalized, use_dimred = "PCA", 
                              n_dimred = 50, ncomponents = 2)

rownames(sce_koh_normalized) <- rownames(sce_koh_normalized)  %>%
  str_replace_all("\\..*", "")

## Bulk RNA data

bulkrna_de <- fread(args$koh_bulkrna, sep = ',') %>%
  dplyr::mutate(geneID=str_replace_all(geneID, "\\..*", ""))

de_mat <- bulkrna_de %>%
  dplyr::select(geneSymbol, geneID, matches("^DE")) %>%
  dplyr::select(-contains("DLL1neg")) %>%
  dplyr::select(-contains("Cardiac"))

lfc_mat <- bulkrna_de %>%
  dplyr::select(geneSymbol, geneID, matches("^sLFC")) %>%
  dplyr::select(-contains("DLL1neg")) %>%
  dplyr::select(-contains("Cardiac"))

l_rel_expr <- plyr::rbind.fill(lapply(1:nrow(lfc_mat), function(i) {
  row <- lfc_mat[i,]
  idx <- which(str_detect(colnames(row), "^sLFC"))
  lfc_cols <- colnames(row)[idx]
  
  types <- strsplit(str_replace(lfc_cols, "^sLFC ", ""), "\\-")
  types_t <- types %>% lapply(
    function(x) return(unname(koh_common_celltype_map[x]))
  )
  
  exprs <- rep(NA, length(koh_bulk_scrna_celltypes))
  names(exprs) <- koh_bulk_scrna_celltypes
  exprs["hESC"] <- 0
  j <- 0
  
  while (any(is.na(exprs))) {
    known_types <- names(which(!is.na(exprs)))
    update_idxs <- sapply(types_t, function(x) (x[2] %in% known_types) & !(x[1] %in% known_types))
    base_vals <- unname(exprs[sapply(types_t[update_idxs], function(x) x[2])])
    update_types <- sapply(types_t[update_idxs], function(x) x[1])
    
    exprs[update_types] <- base_vals + as.numeric(unname(as.data.frame(row)[,idx[update_idxs]]))
  }
  
  exprs <- exprs - min(exprs)
  
  data.frame(geneSymbol=row$geneSymbol, geneID=row$geneID, exprs %>% t)
}))

expr_mat <- l_rel_expr %>% 
  dplyr::filter(geneID %in% rownames(sce_koh_normalized)) %>%
  tibble::column_to_rownames(var = "geneID") %>%
  dplyr::select(-c(geneSymbol))

maxdiffs <- apply(expr_mat, 1, function(x) max(diff(sort(x))))
thres_vals <- apply(expr_mat, 1, function(x) sort(x)[which.max(diff(sort(x)))])

expr_mat_thres <- plyr::rbind.fill(lapply(1:nrow(expr_mat), function(i) {
  expr_mat[i,] <- binarize(expr_mat[i,], thres_vals[i])
}))
rownames(expr_mat_thres) <- rownames(expr_mat)

quantile_cutoff <- .8

rho <- expr_mat_thres[maxdiffs >= quantile(maxdiffs, quantile_cutoff),] %>%
  as.matrix

rho <- rho[intersect(rownames(rho), rownames(sce_koh_normalized)),]
B <- 20
s <- sizeFactors(koh_normalized)

res <- cellassign_em(exprs_obj = sce_koh_normalized[rownames(rho),], 
                     s = s, rho = rho, X = NULL, B = B, 
                     use_priors = TRUE, prior_type = "shrinkage", 
                     delta_variance_prior = FALSE, verbose = FALSE, 
                     em_convergence_thres = 1e-5, num_runs = 3, 
                     min_delta = quantile(maxdiffs, quantile_cutoff))


sce_koh_normalized$cellassign_cluster <- res$cell_type

plotUMAP(sce_koh_normalized, colour_by = "cellassign_cluster")

evaluation_measures <- compute_evaluation_measures(sce_koh_normalized, 
                                                   truth_labels = sce_koh_normalized$celltype,
                                                   inferred_labels = sce_koh_normalized$cellassign_cluster)


sce_koh_normalized@metadata$evaluation_measures <- data.frame(clustering_method='cellassign-shrinkage',
                                                              evaluation_measures)
