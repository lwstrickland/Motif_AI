#!/usr/bin/env Rscript
#
# AddDNAMethylation.R
#
# Calculates mean methylation percentages for motifs in a context-aware manner,
# adds percentages as features.
#
# Usage:
#   Rscript AddDNAMethylation.R \
#     --input /path/to/output_dir/feature_table.rds \
#     --CpG CpG_context.bismark.cov.gz \
#     --CHG CHG_context.bismark.cov.gz \
#     --CHH CHH_context.bismark.cov.gz \
#     --mito_id "ChrM" \
#     --plastid_id "ChrC"
#     --chr_col chromosome --start_col start --end_col stop \
#     --output /path/to/output_dir/feature_table.rds

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(glue)
  library(GenomicRanges)
})

# ---- CLI arguments ----------------------------------------------------
p <- ArgumentParser(description = "Add a bigWig signal summary (default statistic: mean) as a feature")

p$add_argument("--input", help = "Path to input feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--CpG", help = "CpG context coverage file", type = "character", required = TRUE)
p$add_argument("--CHG", help = "CHG context coverage file", type = "character", required = TRUE)
p$add_argument("--CHH", help = "CHH context coverage file", type = "character", required = TRUE)
p$add_argument("--mito_id", help = "Mitochondrial genome sequence ID; Cs on this sequence will be filtered out", type = "character", default = "ChrM")
p$add_argument("--plastid_id", help = "Plastid genome sequence ID; Cs on this sequence will be filtered out", type = "character", default = "ChrC")
p$add_argument("--chr_col", help = "Name of the chromosome column in feature_table", type = "character", default = "chromosome")
p$add_argument("--start_col", help = "Name of the start column in feature_table", type = "character", default = "start")
p$add_argument("--end_col", help = "Name of the end column in feature_table", type = "character", default = "stop")
p$add_argument("--output", help = "Path to write the updated feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--model_zero_based", help = "Is feature_table's start column 0-based (BED-style)?", type = "logical", default = TRUE)

argv <- p$parse_args()

# ---- Load data ----------------------------------------------------------
feat_table <- readRDS(argv$input)

message("Reading CpG methylation coverage from: ", argv$CpG)
CpG <- read_tsv(file = argv$CpG,
                col_names = c("chromosome", "start", "stop", "methylation_percentage", "count_methylated", "count_unmethylated"),
                show_col_types = F)

message("Reading CHG methylation coverage from: ", argv$CHG)
CHG <- read_tsv(file = argv$CHG,
                col_names = c("chromosome", "start", "stop", "methylation_percentage", "count_methylated", "count_unmethylated"),
                show_col_types = F)

message("Reading CHH methylation coverage from: ", argv$CHH)
CHH <- read_tsv(file = argv$CHH,
                col_names = c("chromosome", "start", "stop", "methylation_percentage", "count_methylated", "count_unmethylated"),
                show_col_types = F)

# ---- Filter -----------------------------------------------------------------
# Build function
filter_methylation <- function(x) {
  x_filter <- x |>
    filter(!str_detect(chromosome, paste0(argv$mito_id, "|", argv$plastid_id))) |> # Remove Cs on mitochondria and chloroplast
    mutate(total_count = count_methylated + count_unmethylated) # make new col for total num. reads accounting for each C
  p99 <- quantile(x_filter$total_count, 0.999) # get total read count corresponding to 99.9th percentile
  x_filter <- x_filter |>
    filter(between(total_count, 5, p99)) # filter out Cs w/ less than 5 total supporting reads & Cs in/above 99.9th percentile for total supporting reads
  return(x_filter)
}

# Call function
CpG_filter <- filter_methylation(x = CpG)
CHG_filter <- filter_methylation(x = CHG)
CHH_filter <- filter_methylation(x = CHH)

# ---- Convert to GRanges -----------------------------------------------------
# Model table: Input is BED-based: +1 to start position for GRanges
model_start_offset <- if (argv$model_zero_based) 1L else 0L

gr_model <- GRanges(
  seqnames = feat_table[[argv$chr_col]],
  ranges = IRanges(start = feat_table[[argv$start_col]] + model_start_offset,
                   end = feat_table[[argv$end_col]]),
  motif = feat_table$motif
)

# Methylation bismark.cov: Inputs are not BED-based: start position remains as is for GRanges
gr_CpG <- GRanges(
  seqnames = CpG_filter$chromosome,
  ranges = IRanges(start = CpG_filter$start, end = CpG_filter$stop),
  methylation_percentage = CpG_filter$methylation_percentage
)

gr_CHG <- GRanges(
  seqnames = CHG_filter$chromosome,
  ranges = IRanges(start = CHG_filter$start, end = CHG_filter$stop),
  methylation_percentage = CHG_filter$methylation_percentage
)

gr_CHH <- GRanges(
  seqnames = CHH_filter$chromosome,
  ranges = IRanges(start = CHH_filter$start, end = CHH_filter$stop),
  methylation_percentage = CHH_filter$methylation_percentage
)

# ---- Overlap, update model table --------------------------------------------
# Write function
overlap_methyl <- function(model, methyl, base_table) {
  methyl_name <- deparse(substitute(methyl))
  context <- paste0(str_extract(methyl_name, "(?<=_).*"))
  message("Finding overlaps with ", context, ".")
  # Perform overlap, store unique motifs overlapping methylation context
  ov <- findOverlaps(model, methyl)
  hits <- unique(queryHits(ov))
  message("Found ", length(hits), " feat_table rows overlapping ", context, " contexts.")
  message(length(model) - length(hits), " in feat_table are NA (not overlapped with ", context, " contexts).")
  
  # For each motif, take mean of methylation percentage for all overlapping Cs, store as "means" (grouped by motif)
  means <- tapply(mcols(methyl)$methylation_percentage[subjectHits(ov)], queryHits(ov), mean, na.rm = TRUE)
  
  # In model table GRanges, make new column
  mcols(model)[[paste0(context, "_MethylationPercentage")]] <- NA_real_
  mcols(model)[[paste0(context, "_MethylationPercentage")]][as.integer(names(means))] <- as.numeric(means)
  
  # Take updated gr_model --> Update base_table with new methylation percentage column
  upd_feat_table <- base_table |>
    left_join(as_tibble(mcols(model)[, c("motif", paste0(context, "_MethylationPercentage"))]),
              by = "motif")
  return(upd_feat_table)
}

# Call function, chaining each result into the next so all three columns accumulate
upd_feat_table <- overlap_methyl(model = gr_model, methyl = gr_CpG, base_table = feat_table) # CpG
upd_feat_table <- overlap_methyl(model = gr_model, methyl = gr_CHG, base_table = upd_feat_table) # CHG
upd_feat_table <- overlap_methyl(model = gr_model, methyl = gr_CHH, base_table = upd_feat_table) # CHH

# ---- Add missingness indicators ---------------------------------------------
# Build function
add_missing_ind <- function(table, context) {
  # create symbol for column to be read
  read_col <- sym(glue("{context}_MethylationPercentage"))
  # create symbol for new column to be written
  write_col <- sym(glue("{context}_Missing"))
  
  # When read_col is NA (missing), mark write_col as 1; otherwise, mark write_col as 0
  table |> 
    mutate(
      # use := (walrus operator) for the new column name (forces function to realize that write_col and read_col are not external strings/vectors, but rather column names in feature table)
      !!write_col := case_when(
        is.na(!!read_col) ~ 1,
        TRUE ~ 0
      )
    ) |>
    # relocate new column
    relocate(!!write_col, .after = !!read_col) |>
    # 
    mutate(
      !!read_col := case_when(
        is.na(!!read_col) ~ 0,
        TRUE ~ !!read_col
      ))
}

# Call function
message("Adding CpG/CHG/CHH missingness indicators . . .")
upd_feat_table <- add_missing_ind(table = upd_feat_table, context = "CpG")
upd_feat_table <- add_missing_ind(table = upd_feat_table, context = "CHG")
upd_feat_table <- add_missing_ind(table = upd_feat_table, context = "CHH")

# ---- Save output ------------------------------------------------------------
saveRDS(upd_feat_table, argv$output)







