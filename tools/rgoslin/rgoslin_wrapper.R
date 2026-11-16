#!/usr/bin/env Rscript

# Parse command line arguments input, column, full and summary output
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    opts <- list(
        input = NULL,
        column = "Lipid_Name",
        output_full = "output_full.tsv",
        output_summary = "output_summary.tsv"
    )

    # Substitute template arguments with command line arguments
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--input" && (i + 1) <= length(args)) {
            opts$input <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--column" && (i + 1) <= length(args)) {
            opts$column <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--output_full" && (i + 1) <= length(args)) {
            opts$output_full <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--output_summary" && (i + 1) <= length(args)) {
            opts$output_summary <- args[i + 1]
            i <- i + 2
        } else {
            i <- i + 1
        }
    }
    return(opts)
}

opts <- parse_args(args)

# Check for mandatory input
if (is.null(opts$input)) {
    stop("Please define the input file.", call. = FALSE)
}

# Read input data
df <- read.delim(opts$input, check.names = FALSE, stringsAsFactors = FALSE)

# Resolve target column
target_col <- opts$column
col_idx <- suppressWarnings(as.integer(target_col))

if (!is.na(col_idx)) {
    if (col_idx > 0 && col_idx <= ncol(df)) {
        target_col <- colnames(df)[col_idx]
    } else {
        stop(paste0("Column index ", col_idx, " is out of bounds."), call. = FALSE)
    }
}

# Check if target column exists in dataset
if (!target_col %in% colnames(df)) {
    if ("--column" %in% args) {
        stop(
            paste0("Specified column '", target_col, "' not found in input dataset."),
            call. = FALSE
        )
    } else {
        stop(
            "Default column 'Lipid_Name' was not found in the input file. Please specify the target column using --column.",
            call. = FALSE
        )
    }
}

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
write.table(full_results, file = opts$output_full, sep = "\t", quote = FALSE, row.names = FALSE)

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
write.table(summary_results, file = opts$output_summary, sep = "\t", quote = FALSE, row.names = FALSE)
