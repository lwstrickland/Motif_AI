#!/usr/bin/env Rscript
#
# GenomicFeatureOverlap.R
#
# Overlaps motifs with user-defined genomic feature ranges in BED format.
#
# Usage:
#   Rscript GenomicFeatureOverlap.R \
#     --input /path/to/output_dir/feature_table.rds \
#     --bed features.bed \
#     --column new_column_name \
#     --chr_col chromosome --start_col start --end_col stop \
#     --output /path/to/output_dir/feature_table.rds

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(GenomicRanges)
})

# ---- CLI arguments -----------------------------------------------------------
p <- ArgumentParser(description = "Add overlaps with genomic feature ranges as features")

p$add_argument("--input", help = "Path to input feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--bed", help = "Path to genomic features BED file", type = "character", required = TRUE)
p$add_argument("--column", help = "Name of new overlap column to create", type = "character", required = TRUE)
p$add_argument("--chr_col", help = "Name of the chromosome column in feature_table", type = "character", default = "chromosome")
p$add_argument("--start_col", help = "Name of the start column in feature_table", type = "character", default = "start")
p$add_argument("--end_col", help = "Name of the end column in feature_table", type = "character", default = "stop")
p$add_argument("--output", help = "Path to write updated feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--model_zero_based", help = "Is feature_table's start column 0-based (BED-style)? If TRUE, 1 is added when building GRanges.", type = "logical", default = TRUE)

argv <- p$parse_args()

# ---- Load data ---------------------------------------------------------------
feat_table <- readRDS(argv$input)

message("Reading overlapping features BED from: ", argv$bed)
other_features <- read_tsv(
  argv$bed,
  col_names = c("chr", "start", "end"),
  col_types = cols(chr = col_character(), start = col_double(), end = col_double(), .default = col_guess()),
  comment = "#"
)
# If `other_features` is lowercase in chr column, convert first letter to uppercase to match feature_table's chromosome naming convention.
if (all(grepl("^[a-z]", other_features$chr))) {
  message("Converting first letter of chromosome names in other_features to uppercase to match feature_table.")
  other_features <- other_features |>
    mutate(chr = sub("^([a-z])", "\\U\\1", chr, perl = TRUE))
}

# ---- Initialize overlap column -----------------------------------------------
if (!argv$column %in% colnames(feat_table)) {
  message("Column '", argv$column, "' not found — initializing to 0.")
  feat_table <- feat_table %>% mutate(!!argv$column := 0L)
} else {
  message("Column '", argv$column, "' already exists — incrementing existing values.")
}

# ---- Build GRanges objects ---------------------------------------------------
# If feature_table's start coordinate is 0-based (BED-derived), shift by +1 to
# get 1-based coordinates for GRanges. The end/stop coordinate is unchanged
# either way, since BED's half-open end equals the 1-based inclusive end.
start_offset <- if (argv$model_zero_based) 1L else 0L

gr_model <- GRanges(
  seqnames = feat_table[[argv$chr_col]],
  ranges = IRanges(start = feat_table[[argv$start_col]] + start_offset,
                   end = feat_table[[argv$end_col]])
)

# BED format is 0-based, half-open -> convert other_features start to 1-based for GRanges.
gr_other <- GRanges(
  seqnames = other_features$chr,
  ranges = IRanges(start = other_features$start + 1L, end = other_features$end)
)

# ---- Find overlaps and increment counts ---------------------------------
message("Finding overlaps...")
hits <- unique(queryHits(findOverlaps(gr_model, gr_other)))
message("Found ", length(hits), " feature_table rows overlapping features.")

feat_table[[argv$column]][hits] <- feat_table[[argv$column]][hits] + 1L

# ---- Save output ----------------------------------------------------------
saveRDS(feat_table, argv$output)






