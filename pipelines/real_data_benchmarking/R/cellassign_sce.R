#' Cell type assignments with CellAssign on a SingleCellExperiment

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
library(cellassign)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Run CellAssign on a SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character', nargs = '+',
                    help="Path(s) to SingleCellExperiment RDS")
parser$add_argument('--marker_gene_matrix', metavar='FILE', type = 'character',
                    help="Path to a rho matrix", default = NULL)
parser$add_argument('--marker_list', metavar='FILE', type = 'character',
                    help="Path to a marker list yaml", default = NULL)
parser$add_argument('--include_other', type = 'character',
                    help="Include an other column in the marker gene matrix if using a marker_list")
parser$add_argument('--num_runs', type = 'double',
                    help="Number of runs", default = 1)
parser$add_argument('--rbf_pieces', type = 'integer',
                    help="Number of pieces to fit (RBF)", default = 20)
parser$add_argument('--min_delta', type = 'double',
                    help="Minimum delta value", default = 2)
parser$add_argument('--shrinkage', action = 'store_true',
                    help="Use delta shrinkage prior")
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--design_formula', type = 'character',
                    help="Design matrix formula, or 'none' for none", default = "none")
parser$add_argument('--celltypes', type = 'character', nargs = '+',
                    help="Celltypes to subset", default = "all")
parser$add_argument('--celltype_col', type = 'character',
                    help="Celltype column", default = "celltype")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for cell cycle assignments.")
args <- parser$parse_args()

sce_paths <- unlist(args$sce)

if (length(sce_paths) == 1) {
  sce <- readRDS(sce_paths[1])
} else {
  sces <- lapply(sce_paths, function(x) {
    sce <- readRDS(x)
    colData(sce)$path <- x
    return(sce)
  })
  
  coldata_union <- Reduce(union, lapply(sces, function(x) {
    colnames(colData(x))
  }))
  
  common_genes <- Reduce(intersect, lapply(sces, function(x) {
    rownames(x)
  }))
  
  sces <- lapply(sces, function(x) {
    missing_cols <- setdiff(coldata_union, colnames(colData(x)))
    rowdat_cols <- setdiff(colnames(rowData(x)), c("mean_counts", "log10_mean_counts",
                                                   "n_cells_by_counts", "pct_dropout_by_counts", "total_counts",
                                                   "log10_total_counts"))
    rowData(x) <- rowData(x)[,rowdat_cols]
    colData(x)[,missing_cols] <- NA
    colData(x) <- colData(x)[,coldata_union]
    x <- x[common_genes,]
    return(x)
  })
  
  sce <- Reduce(cbind, sces)
  
  # Recompute size factors
  qclust <- quickCluster(sce, min.size = 30)
  sce <- computeSumFactors(sce, clusters = qclust)
  
  sce$size_factor <- sizeFactors(sce)
}


if (all(unlist(args$celltypes) != "all")) {
  sce$celltype <- as.data.frame(colData(sce))[,args$celltype_col]
  
  sce <- sce %>%
    scater::filter(celltype %in% unlist(args$celltypes))
  
  print(sce)
}

# Attempt to snakemake's default
Sys.setenv(PYTHONPATH='')

reticulate::use_condaenv(args$conda_env, conda = "/home/rstudio/miniconda/bin/conda")

# Process marker gene matrix
if (is.null(args$marker_gene_matrix)) {
  if (is.null(args$marker_list)) {
    stop("Must specify either a marker gene matrix or a marker list file.")
  }
  
  marker_list <- read_yaml(args$marker_list)
  rho <- marker_list_to_mat(marker_list, include_other = as.logical(args$include_other))
} else {
  rho <- read.table(args$marker_gene_matrix, 
                    row.names = 1, 
                    header = TRUE, 
                    sep = "\t",
                    stringsAsFactors = FALSE)
}

rho <- as.matrix(rho)
rho <- rho[intersect(rownames(rho), rowData(sce)$Symbol),]
s <- sizeFactors(sce)

sce_markers <- sce[get_ensembl_id(rownames(rho), sce),]
rownames(sce_markers) <- rownames(rho)
counts(sce_markers) <- as.matrix(counts(sce_markers))

# Create design matrix
if (args$design_formula != "none") {
  dat <- colData(sce_markers) %>% 
    as.data.frame
  design <- model.matrix(as.formula(args$design_formula), data = dat)
} else {
  design <- NULL
}

print(head(design))

cellassign_res <- cellassign(exprs_obj = sce_markers, 
                             s = s, 
                             marker_gene_info = rho, 
                             X = design,
                             B = args$rbf_pieces, 
                             shrinkage = args$shrinkage,  
                             verbose = FALSE, 
                             rel_tol_em = 1e-5, 
                             num_runs = args$num_runs, 
                             min_delta = args$min_delta)

# Write outputs
saveRDS(cellassign_res, file = args$outfname)

cat("Completed.\n")


