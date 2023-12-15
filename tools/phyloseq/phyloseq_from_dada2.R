#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list(
    make_option(c("--sequence_table"), action = "store", dest = "sequence_table", help = "Input sequence table"),
    make_option(c("--taxonomy_table"), action = "store", dest = "taxonomy_table", help = "Input taxonomy table"),
    make_option(c("--output"), action = "store", dest = "output", help = "RDS output")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
# The input sequence_table is an integer matrix
# stored as tabular (rows = samples, columns = ASVs).
seq_table_numeric_matrix <- data.matrix(read.table(opt$sequence_table, sep = "\t"))
# The input taxonomy_table is a table containing
# the assigned taxonomies exceeding the minBoot
# level of bootstrapping confidence. Rows correspond
# to sequences, columns to taxonomic levels. NA
# indicates that the sequence was not consistently
# classified at that level at the minBoot threshold.
tax_table_matrix <- as.matrix(read.table(opt$taxonomy_table, header = FALSE, sep = "\t"))
# Construct a tax_table object.  The rownames of
# tax_tab must match the OTU names (taxa_names)
# of the otu_table defined below.
tax_tab <- tax_table(tax_table_matrix)
# Construct an otu_table object.
otu_tab <- otu_table(seq_table_numeric_matrix, taxa_are_rows = TRUE)
# Construct a phyloseq object.
phyloseq_obj <- phyloseq(otu_tab, tax_tab)
saveRDS(phyloseq_obj, file = opt$output, compress = TRUE)
