#!/usr/bin/env Rscript
#
# GeneRelativePositions.R
#
# Adds motif's distance to nearest gene and position relative to gene as feature.
#
# Usage:
#   Rscript GeneRelativePositions.R \
#     --input /path/to/output_dir/feature_table.rds \
#     --gff annotation.gff3 \
#     --distance_col Distance_NearestGene \
#     --promoter_dist 1000 --downstream_dist 1000 \
#     --chr_col chromosome --start_col start --end_col stop \
#     --output /path/to/output_dir/feature_table.rds

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
})

# ---- CLI arguments ----------------------------------------------------
p <- ArgumentParser(description = "Adds distance-to-nearest-gene (bp) and gene-relative-region overlap columns to feat_table, using a GFF3 genome annotation")

p$add_argument("--input", help = "Path to model feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--gff", help = "Path to GFF3 genome annotation file (may be gzipped)", type = "character", required = TRUE)
p$add_argument("--distance_col", help = "Name of the new distance-to-nearest-gene column (bp)", type = "character", default = "Distance_NearestGene")
p$add_argument("--promoter_dist", help = "Promoter size (bp) upstream of the transcription start site", type = "integer", default = 1000L)
p$add_argument("--downstream_dist", help = "Downstream region size (bp) past the transcription termination site", type = "integer", default = 1000L)
p$add_argument("--gene_feature_types", help = "Comma-separated GFF3 type(s) used to define genes (for distance/promoter/downstream)", type = "character", default = "gene")
p$add_argument("--chr_col", help = "Name of the chromosome column", type = "character", default = "chromosome")
p$add_argument("--start_col", help = "Name of the start column", type = "character", default = "start")
p$add_argument("--end_col", help = "Name of the end column", type = "character", default = "stop")
p$add_argument("--output", help = "Path to write the updated feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--model_zero_based", help = "Is feat_table's start column 0-based (BED-style)? If TRUE, 1 is added when building GRanges.", type = "logical", default = TRUE)

argv <- p$parse_args()

# ---- Load data --------------------------------------------------------------
feat_table <- readRDS(argv$input)

message("Reading GFF3 annotation from: ", argv$gff)
gff <- import(argv$gff, format = "gff3")
# Clean up gff GRanges object (remove unnecessary/uninformative columns)
mcols(gff) <- mcols(gff)[, c("type", "ID", "locus_type", "Parent")]
# Remove rows of genes on mitochondrial/plastid genomes
gff <- gff[!seqnames(gff) %in% c("ChrM", "ChrC")]

# ---- Build feat_table GRanges ---------------------------------------------
model_start_offset <- if (argv$model_zero_based) 1L else 0L

gr_model <- GRanges(
  seqnames = feat_table[[argv$chr_col]],
  ranges = IRanges(start = feat_table[[argv$start_col]] + model_start_offset,
                   end = feat_table[[argv$end_col]])
)

# ---- Build gene-relative region GRanges -------------------------------------
gene_types <- str_split(argv$gene_feature_types, ",")[[1]] |> str_trim() # for if you want to utilize multiple gene feature types (e.g., "gene" and "mRNA"; default is just "gene")
message("Using GFF3 type(s) [", paste(gene_types, collapse = ", "), "] to define genes.")
gr_gene <- gff[gff$type %in% gene_types]
if (length(gr_gene) == 0) {
  stop("No features of type(s) '", argv$gene_feature_types, "' found in GFF3.")
}

# 5'-UTR, 3'-UTR, and exon: Regions defined inherently by GFF3 annotation
gr_5utr <- gff[gff$type == "five_prime_UTR"]
gr_3utr <- gff[gff$type == "three_prime_UTR"]
gr_exon <- gff[gff$type == "exon"]

# Promoter: --promoter_dist bp upstream of the TSS. Downstream: --downstream_dist
# bp past the TES. flank() is strand-aware, so this works for both + and -
# strand genes without special-casing.
gr_promoter <- flank(gr_gene, width = argv$promoter_dist, start = TRUE)
start(gr_promoter) <- pmax(1L, start(gr_promoter)) # ensures no region coordinates are negative by taking the larger value between 1 and (TSS_position - 1000)

gr_downstream <- flank(gr_gene, width = argv$downstream_dist, start = FALSE)
start(gr_downstream) <- pmax(1L, start(gr_downstream)) # ensures no region coordinates are negative by taking the larger value between 1 and (TTS_position + 1000)

# Introns: not an explicit GFF3 feature type, so derive them as the gaps
# between each transcript's exons (grouped by the exon's Parent transcript ID).
transcript_id <- vapply(gr_exon$Parent, function(x) if (length(x)) x[1] else NA_character_, character(1)) # for each exon, pull out corresponding transcript ID
exons_by_tx <- split(gr_exon, transcript_id) # create GRangesList; each element in list is GRanges object of exons corresponding to one transcript_id
tx_span <- unlist(range(exons_by_tx)) # range() collapses exons into a single region for each gene (start codon to stop codon); unlist() flattens the GRangesList into a single Granges object
gr_intron <- unlist(psetdiff(tx_span, exons_by_tx)) # psetdiff(): extracts regions in tx_span not currently in exons_by_tx (in other words, it pulls out the intron coordinates), delivering a GRangesList object; unlist() flattens the GRangesList into a single Granges object of introns

# ---- Distance to nearest gene ------------------------------------------------
message("Computing distance to nearest gene...")
d2n <- distanceToNearest(gr_model, gr_gene) # get distance of each motif to nearest gene
feat_table[[argv$distance_col]] <- NA_integer_ # initialize new col in feat_table, fill w/ NA_integer_
feat_table[[argv$distance_col]][queryHits(d2n)] <- mcols(d2n)$distance # for motifs (queryHits(d2n)), get distance to nearest gene (mcols(d2n)$distance) and insert that distance in bp into feat_table
n_na <- sum(is.na(feat_table[[argv$distance_col]])) # do any motifs still report NA_integer_ as distance to gene?
if (n_na > 0) {
  message(n_na, " feat_table rows have no gene on the same sequence and were left NA in '", argv$distance_col, "'.")
}

# ---- Gene-relative region overlap categories --------------------------------
# Build function for getting gene-relative positional overlaps for motifs
overlap_flag <- function(model, features) {
  as.integer(countOverlaps(model, features) > 0)
}

# Call function --> assign overlaps in gene-relative regions (or lack thereof) as new columns in feat_table
message("Finding overlaps with gene-relative regions...")
feat_table$Overlaps_Promoter <- overlap_flag(gr_model, gr_promoter)
feat_table$Overlaps_FivePrimeUTR <- overlap_flag(gr_model, gr_5utr)
feat_table$Overlaps_Intron <- overlap_flag(gr_model, gr_intron)
feat_table$Overlaps_ThreePrimeUTR <- overlap_flag(gr_model, gr_3utr)
feat_table$Overlaps_Downstream <- overlap_flag(gr_model, gr_downstream)
feat_table$Overlaps_Exon <- overlap_flag(gr_model, gr_exon)

for (col in c("Overlaps_Promoter", "Overlaps_FivePrimeUTR", "Overlaps_Intron",
              "Overlaps_ThreePrimeUTR", "Overlaps_Downstream", "Overlaps_Exon"
              )) {
  message("  ", col, ": ", sum(feat_table[[col]]), " / ", nrow(feat_table), " rows")
}

# ---- Save output ----------------------------------------------------------
saveRDS(feat_table, argv$output)




