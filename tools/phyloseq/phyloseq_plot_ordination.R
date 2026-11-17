#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))

option_list <- list(
    make_option(c("--input"), action = "store", dest = "input", help = "Input file containing a phyloseq object"),
    make_option(c("--method"), action = "store", dest = "method", help = "Ordination method"),
    make_option(c("--distance"), action = "store", dest = "distance", help = "Distance method"),
    make_option(c("--type"), action = "store", dest = "type", help = "Plot type"),
    make_option(c("--output"), action = "store", dest = "output", help = "Output")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options
# Construct a phyloseq object.
phyloseq_obj <- readRDS(opt$input)
# Transform data to proportions as appropriate for
# Bray-Curtis distances.
proportions_obj <- transform_sample_counts(phyloseq_obj, function(otu) otu / sum(otu))
ordination_obj <- ordinate(proportions_obj, method = opt$method, distance = opt$distance)
# Start PDF device driver and generate the plot.
dev.new()
pdf(file = opt$output)
plot_ordination(proportions_obj, ordination_obj, type = opt$type)
dev.off()
