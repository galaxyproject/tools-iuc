#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))

option_list <- list(
    make_option(c("--input"), action = "store", dest = "input", help = "Input RDS file containing a phyloseq object"),
    make_option(c("--output"), action = "store", dest = "output", help = "Output PDF")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
phyloseq_obj <- readRDS(opt$input)
# Start PDF device driver and generate the plot.
dev.new()
pdf(file = opt$output)
plot_richness(phyloseq_obj, x = "samples", color = "samples")
dev.off()
