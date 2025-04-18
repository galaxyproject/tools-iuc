suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))

# Define command-line options
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the phyloseq RDS file", metavar = "FILE"),
    make_option(c("-r", "--rank"), type = "character", help = "Taxonomic rank for aggregation"),
    make_option("--exclude_otu_ids", action = "store_true", default = FALSE, help = "Exclude OTU IDs from output"),
    make_option("--single_rank", action = "store_true", default = FALSE, help = "Only output the specified rank column"),
    make_option("--exclude_na_values", action = "store_true", default = FALSE, help = "Exclude NA values during tax_glom")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate arguments
if (is.null(opt$input) || is.null(opt$rank)) {
    stop("Error: --input and --rank are required arguments.")
}

if (opt$single_rank && !opt$exclude_otu_ids) {
    stop("Error: --single_rank can only be used if --exclude_otu_ids is also specified.")
}

# Load the phyloseq object
physeq <- readRDS(opt$input)

# Print available taxonomic ranks
cat("Available taxonomic ranks:\n")
print(rank_names(physeq))

# Print original number of OTUs
cat("Original number of OTUs:", ntaxa(physeq), "\n")

# Perform tax_glom
physeq_agg <- tax_glom(physeq, taxrank = opt$rank, NArm = opt$exclude_na_values)

# Print new number of taxa after agglomeration
cat("Number of taxa after agglomeration at", opt$rank, "level:", ntaxa(physeq_agg), "\n")

# Extract the taxonomy table after agglomeration
tax_table_agg <- as.data.frame(tax_table(physeq_agg))

# Convert taxonomic columns to character to preserve NA values
tax_table_agg[] <- lapply(tax_table_agg, as.character)

# Add OTU ID column unless excluded
if (!opt$exclude_otu_ids) {
    tax_table_agg <- cbind("OTU ID" = rownames(tax_table_agg), tax_table_agg)
}

# Extract OTU abundance table and convert to data frame
otu_table_agg <- as.data.frame(otu_table(physeq_agg))

# Append taxonomic information to output
otu_table_agg <- cbind(tax_table_agg, otu_table_agg)

tax_table_agg <- otu_table_agg

if (opt$single_rank) {
    # Keep only the specified taxonomic rank column and numeric count columns
    tax_table_agg <- tax_table_agg %>% select(all_of(opt$rank), where(is.numeric))

    # Group by taxonomic rank and sum the counts
    tax_table_agg <- tax_table_agg %>%
        group_by(across(all_of(opt$rank))) %>%
        summarise(across(where(is.numeric), sum), .groups = "drop")
}

# Save the output as a TSV file
output_file <- paste0("physeq_", opt$rank, "_table.tsv")
write.table(tax_table_agg, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
