#!/usr/bin/env Rscript

# Parse command line arguments input, column, full and summary output
args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
column <- args[2]
output_full <- args[3]
output_summary <- args[4]

# Read input data
df <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)

# Resolve target column
col_idx <- as.integer(column)

if (is.na(col_idx) || col_idx < 1 || col_idx > ncol(df)) {
  stop("Invalid lipid name column selection.", call. = FALSE)
}

target_col <- colnames(df)[col_idx]

# Load rgoslin
suppressPackageStartupMessages({
    library(rgoslin)
})

lipid_names <- as.character(df[[target_col]])
lipid_names[is.na(lipid_names)] <- ""

parsed_results <- rgoslin::parseLipidNames(
    lipidNames = lipid_names
)

# Improve error messages for users
empty_idx <- is.na(df[[target_col]]) |
    trimws(as.character(df[[target_col]])) == ""

unparsed_idx <- parsed_results$Grammar == "NOT_PARSEABLE" |
    is.na(parsed_results$Normalized.Name)

# Empty input
parsed_results$Message[empty_idx] <- "No lipid name provided."

# Unparseable lipid names
parsed_results$Message[!empty_idx & unparsed_idx] <-
    "Unsupported lipid nomenclature."

# Full output with all generated columns
full_results <- cbind(df, parsed_results)
write.table(full_results, file = output_full, sep = "\t", quote = FALSE, row.names = FALSE)

# Summary output with selected columns
summary_cols <- c(
    "Normalized.Name", "Grammar", "Message",
    "Lipid.Maps.Category", "Lipid.Maps.Main.Class",
    "Species.Name", "Molecular.Species.Name",
    "Level", "Total.C", "Total.DB", "Mass", "Sum.Formula"
)
keep_cols <- intersect(summary_cols, colnames(parsed_results))
if (length(keep_cols) > 0) {
    parsed_summary <- parsed_results[, keep_cols, drop = FALSE]
} else {
    parsed_summary <- parsed_results
}

summary_results <- cbind(df, parsed_summary)
write.table(summary_results, file = output_summary, sep = "\t", quote = FALSE, row.names = FALSE)
