#' Find immune-associated SNVs, indels, and CNVs

library(tidyverse)
library(data.table)
library(methods)
library(cowplot)
library(xlsx)
library(ImportExport)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

library(scrna.utils)
library(scrna.sceutils)
library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Find immune-associated variants by immune term.")
parser$add_argument('--snvs', type = 'character', metavar='FILE',
                    help = "SNV table.")
parser$add_argument('--indels', type = 'character', metavar='FILE',
                    help = "Indel table.")
parser$add_argument('--cnvs', type = 'character', metavar='FILE',
                    help = "CNV table.")
parser$add_argument('--pancancer_annotations', type = 'character', metavar='FILE',
                    help = "Nanostring pancancer panel annotations.")
parser$add_argument('--sample_names', type = 'character', nargs = '+',
                    help = "Samples to consider.")
parser$add_argument('--museq_probability', type = 'double',
                    help = "MutationSeq probability to filter at", default = 0.9)
parser$add_argument('--outdir', type = 'character', metavar = 'DIR',
                    help="Output directory for filtered results")
args <- parser$parse_args()

snvs <- fread(args$snvs)
indels <- fread(args$indels)
cnvs <- fread(args$cnvs)
sample_names <- unlist(args$sample_names)

pancancer_annotations <- fread(args$pancancer_annotations)
immune_gene_table <- pancancer_annotations %>% 
  dplyr::filter(!gene_class %in% c("HK", "CT Antigen"))
additional_genes <- c("CREBBP")
interesting_gene_table <- immune_gene_table %>% 
  plyr::rbind.fill(data.frame(Gene=additional_genes, gene_class = "interesting"))

assembly_gene_list <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

snvs_filtered <- snvs %>%
  dplyr::mutate(sample=paste(patientName, timepoint, sep = "")) %>%
  dplyr::filter(sample %in% sample_names,
                !is.na(geneName),
                PR >= args$museq_probability) %>%
  dplyr::rename(Gene=geneName) %>%
  dplyr::inner_join(interesting_gene_table, by = "Gene") %>%
  dplyr::filter(impact != "MODIFIER")

indels_filtered <- indels %>%
  dplyr::mutate(sample=str_extract(sampleID, "^.*(?=_)")) %>%
  dplyr::filter(sample %in% sample_names,
                IMPACT != "MODIFIER") %>%
  dplyr::rename(Gene=GENE) %>% 
  dplyr::inner_join(interesting_gene_table, by = "Gene")

cnvs_formatted <- cnvs %>%
  dplyr::mutate(sample=str_extract(Sample, "^.*(?=_)"),
                Chromosome = paste0("chr", Chromosome)) %>%
  dplyr::filter(sample %in% sample_names,
                !(MajorCN == 1 & MinorCN == 1)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           ignore.strand = TRUE,
                           seqnames.field = c("Chromosome"),
                           start.field = "Start_Position",
                           end.field = "End_Position")

assembly_gene_list$Symbol <- mapIds(org.Hs.eg.db, assembly_gene_list$gene_id, 'SYMBOL', 'ENTREZID')
assembly_gene_list_interesting <- assembly_gene_list[assembly_gene_list$Symbol %in% interesting_gene_table$Gene]

cnvs_formatted_gene_annotated <- findOverlaps(cnvs_formatted, assembly_gene_list_interesting, ignore.strand = TRUE) %>%
  as.data.frame
cnvs_formatted_gene_annotated$Gene <- assembly_gene_list_interesting$Symbol[cnvs_formatted_gene_annotated$subjectHits]

cnvs_filtered <- cnvs_formatted %>%
  as.data.frame %>%
  dplyr::mutate(queryHits=1:n()) %>%
  dplyr::inner_join(cnvs_formatted_gene_annotated, by = "queryHits") %>%
  dplyr::inner_join(interesting_gene_table, by = "Gene")

# Write outputs

if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

filtered_snv_file <- file.path(args$outdir, "snvs_filtered.tsv")
filtered_indel_file <- file.path(args$outdir, "indels_filtered.tsv")
filtered_cnv_file <- file.path(args$outdir, "cnvs_filtered.tsv")

write.table(snvs_filtered, file = filtered_snv_file,
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)
write.table(indels_filtered, file = filtered_indel_file,
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)
write.table(cnvs_filtered, file = filtered_cnv_file,
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)

cat("Completed.\n")

