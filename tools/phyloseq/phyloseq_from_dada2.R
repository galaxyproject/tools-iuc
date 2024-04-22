#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list(
    make_option(c("--sequence_table"), action = "store", dest = "sequence_table", help = "Input sequence table"),
    make_option(c("--taxonomy_table"), action = "store", dest = "taxonomy_table", help = "Input taxonomy table"),
    make_option(c("--sample_table"), action = "store", default = NULL, dest = "sample_table", help = "Input sample table"),
    make_option(c("--output"), action = "store", dest = "output", help = "RDS output")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
# The input sequence_table is an integer matrix
# stored as tabular (rows = samples, columns = ASVs).
seq_table_numeric_matrix <- data.matrix(read.table(opt$sequence_table, header = T, sep = "\t", row.names = 1, check.names = FALSE))
# The input taxonomy_table is a table containing
# the assigned taxonomies exceeding the minBoot
# level of bootstrapping confidence. Rows correspond
# to sequences, columns to taxonomic levels. NA
# indicates that the sequence was not consistently
# classified at that level at the minBoot threshold.
tax_table_matrix <- as.matrix(read.table(opt$taxonomy_table, header = T, sep = "\t", row.names = 1, check.names = FALSE))
# Construct a tax_table object.  The rownames of
# tax_tab must match the OTU names (taxa_names)
# of the otu_table defined below.
tax_tab <- tax_table(tax_table_matrix)

# Construct an otu_table object.
otu_tab <- otu_table(seq_table_numeric_matrix, taxa_are_rows = TRUE)

# Construct a phyloseq object.
phyloseq_obj <- phyloseq(otu_tab, tax_tab)
if (!is.null(opt$sample_table)) {
    sample_tab <- sample_data(
        read.table(opt$sample_table, header = T, sep = "\t", row.names = 1, check.names = FALSE)
    )
    phyloseq_obj <- merge_phyloseq(phyloseq_obj, sample_tab)
}

# use short names for our ASVs and save the ASV sequences
# refseq slot of the phyloseq object as described in
# https://benjjneb.github.io/dada2/tutorial.html
dna <- Biostrings::DNAStringSet(taxa_names(phyloseq_obj))
names(dna) <- taxa_names(phyloseq_obj)
phyloseq_obj <- merge_phyloseq(phyloseq_obj, dna)
taxa_names(phyloseq_obj) <- paste0("ASV", seq(ntaxa(phyloseq_obj)))

print(phyloseq_obj)

# save R object to file
saveRDS(phyloseq_obj, file = opt$output, compress = TRUE)
