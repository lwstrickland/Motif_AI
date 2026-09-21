#!/usr/bin/env Rscript
#
# ArrangeFeatureTable.R
#
# Arranges initial MotifAI feature table from FIMO (v5.5.4) output.
#
# Usage:
#   Rscript ArrangeFeatureTable.R \
#     --input /path/to/output_dir/FIMO_results/fimo.tsv \
#     --output /path/to/output_dir/feature_table.rds
#     --metadata /path/to/Zenodo_files/OMalley2016_DAPseq_metadata/OMalley2016_DAP_At_genexmotif_FRiPThresh10perc.tsv

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})

# ---- CLI arguments -----------------------------------------------------------

p <- ArgumentParser(description = "Initialize feature table from FIMO output")

p$add_argument("--input", help = "Path to input FIMO_results/fimo.tsv file", type = "character", required = TRUE)
p$add_argument("--output", help = "Path to write feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--metadata", help = "Path to O'Malley 2016 TF/motif metadata table", type = "character", required = TRUE)

argv <- p$parse_args()

# ---- Files -------------------------------------------------------------------
# FIMO results
message("Reading FIMO results from: ", argv$input)
fimo_res <- suppressWarnings(read_tsv(file = argv$input,
                     col_names = TRUE,
                     show_col_types = FALSE))
# Remove FIMO run descriptors at the bottom +
# create new col `motif_base_id`
fimo_res <- fimo_res |>
  slice(1:(n() - 3)) |>
  mutate(motif_base_id = str_remove(motif_id, "\\.\\d+$"))
message("Num. FIMO hits: ", nrow(fimo_res))


# Metadata for motif hits/TF gene IDs for O'Malley 2016 DAP-seq TFs w/ FRiP score of at least 10%
metadata_DAP <- suppressWarnings(read_tsv(file = argv$metadata,
                         col_names = TRUE,
                         show_col_types = FALSE))

# ---- Operations --------------------------------------------------------------
# Join fimo_res w/ metadata_DAP to get some additional helpful metadata +
#  arrange into initial feature table
feature_table <- fimo_res |>
  left_join(metadata_DAP, join_by("motif_base_id" == "base_id")) |> # join
  select(motif_id, motif_alt_id, Gene_ID, class, family, type, sequence_name, start, stop, matched_sequence, score, `p-value`, `q-value`) |> # select desired columns
  rename(motif = motif_id,
         name = motif_alt_id,
         source = type,
         chromosome = sequence_name) |> # rename some cols
  mutate(motif = paste0(motif, "_", chromosome, ":", start, "-", stop), # expand entries in motif, make `motif` names more descriptive
         name = str_remove(name, "^([^.]*\\.){2}")) |> # remove motif ID, keep only corresponding TF name for `name` column
  distinct(motif, .keep_all = TRUE) # remove redundancy in table
message("Num. non-redundant motifs retained: ", nrow(feature_table))

# Convert p- & q-vals into -log10(p) & -log10(q)
message("Converting p- & q-vals into -log10(p) & -log10(q) . . .")
feature_table <- feature_table |>
  mutate(across(.cols = c(`p-value`, `q-value`), ~ -log10(.x), .names = "log10{.col}")) |> # -log10-transform p-values & q-values, create new cols
  relocate(c(`log10p-value`, `log10q-value`), .after = "score") |> # move new cols to after score
  rename("log10p" = "log10p-value",
         "log10q" = "log10q-value") |> # rename
  select(-c(`p-value`, `q-value`)) # remove old p- & q-val columns

# ---- Save feature table ------------------------------------------------------
saveRDS(feature_table, argv$output)





