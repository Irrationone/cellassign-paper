#' Find immune-associated SNVs, indels, and CNVs

library(tidyverse)
library(data.table)
library(methods)
library(cowplot)
library(xlsx)
library(ImportExport)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Find immune-associated variants by immune term.")
parser$add_argument('--pancancer_annotations', type = 'character', metavar='FILE',
                    help = "Nanostring pancancer panel annotations.")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for table")
args <- parser$parse_args()