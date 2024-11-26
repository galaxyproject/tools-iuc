#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("ggplot2"))

# Define options
option_list <- list(
  make_option(c("--input"),
    action = "store", dest = "input",
    help = "Input file containing a phyloseq object"
  ),
  make_option(c("--x"),
    action = "store", dest = "x",
    help = "Variable for x-axis (e.g., 'Sample', 'Phylum')"
  ),
  make_option(c("--fill"),
    action = "store", dest = "fill",
    help = "Variable for fill color (e.g., 'Genus', 'Order')"
  ),
  make_option(c("--facet"),
    action = "store", dest = "facet", default = NULL,
    help = "Facet by variable (optional)"
  ),
  make_option(c("--output"),
    action = "store", dest = "output",
    help = "Output file (PDF)"
  )
)

# Parse arguments
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

# Load phyloseq object
print(paste("Trying to read:", opt$input))
physeq <- readRDS(opt$input)

# Check if the 'x' and 'fill' variables are valid
if (!opt$x %in% colnames(sample_data(physeq))) {
  stop(paste("Error: x variable", opt$x, "does not exist in the sample data."))
}

if (!opt$fill %in% colnames(sample_data(physeq))) {
  stop(paste("Error: fill variable", opt$fill, "does not exist in the sample data."))
}

# Generate bar plot
p <- plot_bar(physeq, x = opt$x, fill = opt$fill)

# Only facet if the facet variable is provided and exists in the sample data
if (!is.null(opt$facet)) {
  if (opt$facet %in% colnames(sample_data(physeq))) {
    p <- p + facet_wrap(as.formula(paste("~", opt$facet)))
  } else {
    warning(paste("Facet variable", opt$facet, "does not exist in the sample data. Faceting will be skipped."))
  }
}

# Save to output file
ggsave(opt$output, plot = p, width = 10, height = 8)
