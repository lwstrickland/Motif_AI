#!/usr/bin/env Rscript
#
# MotifAIPredict.R
#
# Using a constructed feature table + the model weights obtained by training
# XGBoost on feature table, make predictions on motifs -->
# deliver MotifAI output TSV file.
#
# Usage:
#   Rscript MotifAIPredict.R \
#     --input  /path/to/output_dir/feature_table.rds \
#     --hf     model_file.ubj \
#     --output /path/to/output_dir/

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(xgboost)
})

# ---- CLI arguments -----------------------------------------------------------
p <- ArgumentParser(description = "Make predictions on motif hits identified in input DNA sequences")

p$add_argument("--input", help = "Path to model feature_table.rds file", type = "character", required = TRUE)
p$add_argument("--hf", help = "Path to model file, downloaded from Huggingface repo", type = "character", required = TRUE)
p$add_argument("--output", help = "Directory to write the final MotifAI output to", type = "character", required = TRUE)

argv <- p$parse_args()

# ---- Load/prep data ----------------------------------------------------------
feat_table <- readRDS(argv$input)

message("Loading model to make predictions on motifs . . .")
model <- xgb.load(argv$hf)

# Select out only features as data frame
feat_table <- df |>
  select(score:Overlaps_Exon)

# Construct matrix, ensure all inputs are numeric!
feat_matrix <- feat_table |> as.matrix()
class(feat_matrix) <- "numeric"

# Set seed
set.seed(123)

# ---- Predict --> Deliver output! ---------------------------------------------
message("Making predictions . . .")
predictions <- predict(object = model,
                       newdata = feat_matrix)

# Assign predicted probabilities back to individual motifs --> Reorder motifs based on "rank" --> Save MotifAI output!
df_out <- df |>
  mutate(probability = round(predictions, digits = 5)) |>
  arrange(desc(probability)) |>
  mutate(rank = row_number()) |>
  select(motif, name, probability, rank)

# Save output files
## Predictions
write_tsv(x = df_out,
          file = paste0(argv$output, "motifai_Results_", format(Sys.Date(), "%m.%d.%Y"), ".tsv"))

## Full constructed feature table in TSV format (for user's easy viewing on CLI, if they wish to)
write_tsv(x = feat_table,
          file = paste0(argv$output, "feature_table.tsv"))
