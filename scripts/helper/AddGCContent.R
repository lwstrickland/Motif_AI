#!/usr/bin/env Rscript
#
# AddGCContent.R
#
# Adds motif GC content as a feature.
#
# Usage:
#   Rscript AddGCContent.R \
#     --input /path/to/output_dir/feature_table.rds \
#     --output /path/to/output_dir/feature_table.rds

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})

# ---- CLI arguments -----------------------------------------------------------

p <- ArgumentParser(description = "Add GC content as a feature")

p$add_argument("--input", help = "Path to input feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--output", help = "Path to write updated feature_table.rds file", type = "character", required = TRUE)

argv <- p$parse_args()

# ---- File(s) -----------------------------------------------------------------
feat_table <- readRDS(file = argv$input)

# ---- Operations --------------------------------------------------------------
# Calculate per-motif GC content --> Add as feature to table
upd_feat_table <- feat_table |>
  mutate(GC_content = str_count(matched_sequence, "G|C") / str_length(matched_sequence))

# ---- Save feature table ------------------------------------------------------
saveRDS(upd_feat_table, argv$output)



