suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("dplyr"))

# Parse command args
args <- commandArgs(trailingOnly = TRUE)

physeq_file <- args[1]
tax_rank <- args[2]
use_counts <- "--counts" %in% args
exclude_otu_ids <- "--exclude_otu_ids" %in% args
single_rank <- "--single_rank" %in% args
exclude_na_values <- "--exclude_na_values" %in% args

if (single_rank && !exclude_otu_ids) {
  stop("Error: --single can only be used if --exclude_otu_ids is also specified.")
}

# Load the phyloseq object
physeq <- readRDS(physeq_file)

# Print available taxonomic ranks
cat("Available taxonomic ranks:\n")
print(rank_names(physeq))

# Peint original number of OTUs
cat("Original number of OTUs:", ntaxa(physeq), "\n")

# Perform tax_glom
physeq_agg <- tax_glom(physeq, taxrank = tax_rank, NArm = exclude_na_values)

# Print new number of taxa after agglomeration
cat("Number of taxa after agglomeration at", tax_rank, "level:", ntaxa(physeq_agg), "\n")

# Extract the taxonomy table after agglomeration
tax_table_agg <- as.data.frame(tax_table(physeq_agg))

# Convert taxonomic columns to character to preserve NA values
tax_table_agg[] <- lapply(tax_table_agg, as.character)

# Add OTU ID column unless excluded
if (!exclude_otu_ids) {
  tax_table_agg <- cbind("OTU ID" = rownames(tax_table_agg), tax_table_agg)
}

if (use_counts) {
  # Extract OTU abundance table and convert to data frame
  otu_table_agg <- as.data.frame(otu_table(physeq_agg))
  
  # Append taxonomic information to output
  otu_table_agg <- cbind(tax_table_agg, otu_table_agg)
  
  tax_table_agg <- otu_table_agg
}

if (single_rank) {
  # Keep only the specified taxonomic rank column
  tax_table_agg <- tax_table_agg %>% select(all_of(tax_rank), everything())
  
  if (use_counts) {
    # Group by taxonomic rank and sum the counts
    tax_table_agg <- tax_table_agg %>%
      group_by(across(all_of(tax_rank))) %>%
      summarise(across(where(is.numeric), sum), .groups = "drop")
  } else {
    # Remove duplicate values, keeping only unique taxonomic rank entries
    tax_table_agg <- distinct(tax_table_agg, .keep_all = TRUE)
  }
}

# Save the output as a TSV file
output_file <- paste0("physeq_", tax_rank, "_table.tsv")
write.table(tax_table_agg, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)