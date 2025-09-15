#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("tidyverse"))

# Option parsing
option_list <- list(
    make_option(c("--input"),
        action = "store", dest = "input",
        help = "Input file containing a phyloseq object"
    ),
    make_option(c("--output"),
        action = "store", dest = "output",
        help = "Output file for the updated phyloseq object"
    ),
    make_option(c("--ranks"),
        action = "store", dest = "ranks",
        help = "Comma-separated list of taxonomy ranks (default: Kingdom,Phylum,Class,Order,Family,Genus,Species)"
    )
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

cat("Input file: ", opt$input, "\n")
cat("Output file: ", opt$output, "\n")
cat("Ranks provided: ", opt$ranks, "\n")

if (is.null(opt$ranks)) {
    opt$ranks <- "Kingdom,Phylum,Class,Order,Family,Genus,Species"
}

# Parse rank names
rank_names <- unlist(strsplit(opt$ranks, ","))

# Load phyloseq object
physeq <- readRDS(opt$input)

# Check if physeq object is loaded successfully
if (is.null(physeq)) {
    stop("Error: Failed to load the Phyloseq object. Check the input file.")
}

cat("Phyloseq object successfully loaded.\n")
cat("Class of loaded object: ", class(physeq), "\n")

# Check the current tax_table
cat("Current tax_table:\n")
print(tax_table(physeq))


# Strict check for taxonomy table and provided ranks
if (ncol(tax_table(physeq)) != length(rank_names)) {
    stop(
        "Error: Number of columns in tax_table does not match the number of provided ranks. ",
        "Please ensure the taxonomy table matches the ranks exactly."
    )
}

# Set column names to the provided ranks
colnames(tax_table(physeq)) <- rank_names

# Confirm the changes
cat("Updated tax_table:\n")
print(tax_table(physeq))

# Save the updated phyloseq object
saveRDS(physeq, file = opt$output, compress = TRUE)
cat("Updated Phyloseq object saved to: ", opt$output, "\n")
