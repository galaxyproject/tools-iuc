#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("ggplot2"))

option_list <- list(
  make_option(c("--input"),
              action = "store", dest = "input",
              help = "Input file containing a phyloseq object"),
  make_option(c("--x"),
              action = "store", dest = "x",
              help = "Variable for x-axis (e.g., 'Sample', 'Phylum')"),
  make_option(c("--fill"),
              action = "store", dest = "fill",
              help = "Variable for fill color (e.g., 'Genus', 'Order')"),
  make_option(c("--facet"),
              action = "store", dest = "facet", default = NULL,
              help = "Facet by variable (optional)"),
  make_option(c("--output"),
              action = "store", dest = "output",
              help = "Output file (PDF)")
)

# Parse arguments
parser <- OptionParser(
  usage = "%prog [options] file",
  option_list = option_list
)
args <- parse_args(parser)
opt <- args$options

# Load phyloseq object
physeq <- readRDS(opt$input)

# Generate bar plot
p <- plot_bar(physeq, x = opt$x, fill = opt$fill)

if (!is.null(opt$facet)) {
  p <- p + facet_wrap(as.formula(paste("~", opt$facet)))
}

# Save to output file
ggsave(opt$output, plot = p, width = 10, height = 8)
