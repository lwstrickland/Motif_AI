#!/usr/bin/env Rscript
#
# AddBigWigSignal.R
#
# Add signal from BigWig tracks, averaged over the length of motifs.
#
# Usage:
#   Rscript AddBigWigSignal.R \
#     --input /path/to/output_dir/feature_table.rds \
#     --bigwig features.bw \
#     --stat mean \
#     --column new_column_name \
#     --chr_col chromosome --start_col start --end_col stop \
#     --output /path/to/output_dir/feature_table.rds

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
})

# ---- CLI arguments ----------------------------------------------------
p <- ArgumentParser(description = "Add a bigWig signal summary (default statistic: mean) as a feature")

p$add_argument("--input", help = "Path to input feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--bigwig", help = "Path to bigWig (.bw) file", type = "character", required = TRUE)
p$add_argument("--column", help = "Name of the new signal column to create", type = "character", default = "bw_signal")
p$add_argument("--stat", help = "Summary statistic to compute over each feature", type = "character",
               default = "mean", choices = c("mean", "min", "max", "sum", "sd", "coverage"))
p$add_argument("--chr_col", help = "Name of the chromosome column in feature_table", type = "character", default = "chromosome")
p$add_argument("--start_col", help = "Name of the start column in feature_table", type = "character", default = "start")
p$add_argument("--end_col", help = "Name of the end column in feature_table", type = "character", default = "stop")
p$add_argument("--output", help = "Path to write the updated feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--model_zero_based", help = "Is feature_table's start column 0-based (BED-style)?", type = "logical", default = TRUE)

argv <- p$parse_args()

# ---- Load data ----------------------------------------------------------
feat_table <- readRDS(argv$input)

# ---- Build GRanges --------------------------------------------------------
model_start_offset <- if (argv$model_zero_based) 1L else 0L

gr_model <- GRanges(
  seqnames = feat_table[[argv$chr_col]],
  ranges = IRanges(start = feat_table[[argv$start_col]] + model_start_offset,
                   end = feat_table[[argv$end_col]])
)
# Carry the original row index as metadata so we can map results back
# correctly even if some ranges get dropped (e.g. missing chromosomes).
mcols(gr_model)$row_id <- seq_len(nrow(feat_table))

# ---- Open bigWig and align seqlevels ------------------------------------
message("Opening bigWig file: ", argv$bigwig)
bw <- BigWigFile(argv$bigwig)
bw_seqinfo <- seqinfo(bw)

# Detect whether the bigWig's actual chromosome names are lowercase (e.g.
# chr1) before any renaming below -- this reflects the naming the bigWig file
# itself uses internally, which we still need later when querying it directly
# via import() (see the chunk loop below).
bw_is_lowercase <- all(grepl("^[a-z]", seqnames(bw_seqinfo)))

# Fix chromosome names, if need be (chr1 --> Chr1, chr2 --> Chr2. etc.)
if (bw_is_lowercase) {
  message("Converting first letter of bigWig chromosome names to uppercase to match feat_table.")
  seqlevels(bw_seqinfo) <- sub("^([a-z])", "\\U\\1", seqlevels(bw_seqinfo), perl = TRUE)
}

# Restrict to chromosomes present in both feat_table and the bigWig, and
# trim any out-of-bound ranges (protects against off-by-one / assembly
# mismatches causing summary() to error out).
common_chroms <- intersect(seqlevels(gr_model), seqnames(bw_seqinfo))
if (length(common_chroms) == 0) {
  stop("No chromosome names in feat_table match those in the bigWig file. ",
       "Check naming conventions (e.g. 'Chr1' vs 'chr1' vs '1').")
}
if (length(common_chroms) < length(seqlevels(gr_model))) {
  missing_chroms <- setdiff(seqlevels(gr_model), common_chroms)
  warning("The following chromosomes in feat_table are absent from the bigWig ",
          "and will get NA signal values: ", paste(missing_chroms, collapse = ", "))
}

seqlevels(gr_model, pruning.mode = "coarse") <- common_chroms
seqinfo(gr_model) <- bw_seqinfo[common_chroms]
gr_model <- trim(gr_model)

# ---- Compute signal summary ------------------------------------------------
message("Computing '", argv$stat, "' signal over ", length(gr_model), " regions...")
# Process motifs in chunks to avoid excessive memory use
chunk_size <- 50000L # process 50,000 motifs at a time

# Create numeric vector of length(gr_model)
signal <- numeric(length(gr_model))
signal[] <- NA_real_

starts <- seq(1L, length(gr_model), by = chunk_size)

for (i in seq_along(starts)) {
  
  s <- starts[i]
  e <- min(s + chunk_size - 1L, length(gr_model)) # determine e (end) as the end of the next chunk, or the end of the motif list, whichever is less (determined by min())
  
  message(sprintf(
    "  Processing regions %d-%d of %d",
    s, e, length(gr_model)
  ))
  
  chunk <- gr_model[s:e]
  
  # CRITICAL: import() queries the bigWig file's index directly using the
  # seqnames of `which`, so they must match the file's actual chromosome
  # naming -- not gr_model's seqinfo (which was only aligned above for
  # bookkeeping/common_chroms purposes). Only convert to lowercase if the
  # bigWig itself used lowercase names (detected earlier, before bw_seqinfo
  # was capitalized to match feat_table).
  if (bw_is_lowercase) {
    seqlevels(chunk) <- sub("^([A-Z])", "\\L\\1", seqlevels(chunk), perl = TRUE)
  }
  
  # Process motifs in chunks, calculating desired bigWig statistic (mean, min, max, etc.) for each interval
  vals <- import(
    bw,
    which = chunk,
    as = "NumericList"
  )
  
  signal[s:e] <- vapply(
    vals,
    function(x) {
      
      if (length(x) == 0)
        return(NA_real_)
      
      switch(
        argv$stat,
        mean = mean(x),
        max = max(x),
        min = min(x),
        sum = sum(x),
        sd = if (length(x) > 1) sd(x) else 0,
        coverage = sum(x != 0),
        stop("Unknown statistic")
      )
      
    },
    numeric(1)
  )
}

stopifnot(length(signal) == length(gr_model))

# Check
message("Signal summary:")
summary(signal)

uniq_sig_vals <- length(unique(signal))
message("Unique signal values: ", uniq_sig_vals)

num_na_vals <- sum(is.na(signal))
message("Number of NA values: ", num_na_vals)

# ---- Map values back onto the full (untrimmed) feat_table ---------------
# Use the row_id metadata carried through gr_model, rather than assuming
# positional alignment -- this stays correct even if seqlevels pruning or
# trim() dropped/reordered any ranges (e.g. missing chromosomes).
feat_table[[argv$column]] <- NA_real_
feat_table[[argv$column]][mcols(gr_model)$row_id] <- signal

# ---- Save output ----------------------------------------------------------
saveRDS(feat_table, argv$output)






