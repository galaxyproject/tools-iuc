#!/usr/bin/env Rscript

# Load rgoslin
suppressPackageStartupMessages({
    library(rgoslin)
})

# Parse command line arguments input, column, full and summary output
args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
column <- as.integer(args[2])
has_header <- as.logical(args[3])
full_output <- as.logical(args[4])
output_full <- args[5]
output_summary <- args[6]

# Read input data
df <- read.delim(
    input,
    header = has_header,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Resolve target column
if (is.na(column) || column < 1 || column > ncol(df)) {
    stop("Invalid lipid name column selection.", call. = FALSE)
}

lipid_names <- as.character(df[[column]])
lipid_names[is.na(lipid_names)] <- ""

parsed_results <- parseLipidNames(
    lipidNames = lipid_names
)

empty_idx <- is.na(df[[column]]) |
    trimws(as.character(df[[column]])) == ""

parsed_results$Message[empty_idx] <- "No lipid name provided."

# Full output with all generated columns
if (full_output) {
    full_results <- cbind(df, parsed_results)
    write.table(
        full_results,
        file = output_full,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
    )
}

# Summary output with selected columns
summary_cols <- c(
    "Normalized.Name", "Grammar", "Message",
    "Lipid.Maps.Category", "Lipid.Maps.Main.Class",
    "Species.Name", "Molecular.Species.Name",
    "Level", "Total.C", "Total.DB", "Mass", "Sum.Formula"
)

parsed_summary <- parsed_results[, summary_cols, drop = FALSE]

summary_results <- cbind(df, parsed_summary)
write.table(
    summary_results,
    file = output_summary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
)
