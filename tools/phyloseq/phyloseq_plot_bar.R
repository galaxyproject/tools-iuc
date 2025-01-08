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
        action = "store", dest = "fill", default = NULL,
        help = "Variable for fill color (e.g., 'Genus', 'Order') (optional)"
    ),
    make_option(c("--facet"),
        action = "store", dest = "facet", default = NULL,
        help = "Facet by variable (optional)"
    ),
    make_option(c("--output"),
        action = "store", dest = "output",
        help = "Output file (PDF)"
    ),
    make_option(c("--topX"),
        action = "store", dest = "topX", default = NULL,
        help = "Show only the top X taxa based on abundance (e.g., '10') (optional)"
    ),
    make_option(c("--keepOthers"),
        action = "store_true", dest = "keepOthers", default = FALSE,
        help = "Keep taxa not in topX and label them as 'Others' (optional)"
    ),
    make_option(c("--normalize"),
        action = "store_true", dest = "normalize", default = FALSE,
        help = "Normalize abundances to sum to 100% (optional)"
    )
)

# Parse arguments
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

# Validate required options
if (is.null(opt$input) || opt$input == "") {
    stop("Error: Input file is required.")
}
if (is.null(opt$output) || opt$output == "") {
    stop("Error: Output file is required.")
}

# Load phyloseq object
print(paste("Trying to read:", opt$input))
physeq <- readRDS(opt$input)

print(rank_names(physeq))
print(opt$fill)

# Check if the 'x' and 'fill' variables are valid
sample_vars <- colnames(sample_data(physeq))

# If topX is provided, filter the phyloseq object to show only top X taxa
if (!is.null(opt$topX) && opt$topX != "") {
    topX <- as.numeric(opt$topX)
    if (is.na(topX) || topX <= 0) {
        stop("Error: topX should be a positive number.")
    }

    # Aggregate the data at the selected rank (e.g., Phylum)
    tax_rank <- opt$fill # Adjust as necessary
    physeq_agg <- tax_glom(physeq, taxrank = tax_rank)

    # Get the abundance of each taxon at the selected rank
    taxa_abundance <- taxa_sums(physeq_agg)

    # Identify the top X taxa
    top_taxa <- names(sort(taxa_abundance, decreasing = TRUE))[1:topX]

    if (opt$keepOthers) {
        # Label taxa not in top_taxa as "Others"
        tax_table(physeq_agg)[, tax_rank][!rownames(tax_table(physeq_agg)) %in% top_taxa] <- "Others"
        physeq <- physeq_agg
    } else {
        # Subset the phyloseq object to keep only the top X taxa
        physeq_filtered <- prune_taxa(top_taxa, physeq_agg)
        physeq <- physeq_filtered
    }
}

# Normalize to relative abundances if requested
if (opt$normalize) {
    print("Normalizing abundances to sum to 100%...")
    physeq <- transform_sample_counts(physeq, function(x) 100 * x / sum(x))
}

# Generate bar plot
if (!is.null(opt$x) && opt$x != "") {
    p <- plot_bar(physeq, x = opt$x, fill = opt$fill)
} else {
    p <- plot_bar(physeq, fill = opt$fill) # If no x is provided, don't include x
}

# Only facet if the facet variable is provided and exists in the sample data
if (!is.null(opt$facet) && opt$facet != "") {
    if (opt$facet %in% sample_vars) {
        p <- p + facet_wrap(as.formula(paste("~", opt$facet)))
    } else {
        warning(paste("Facet variable", opt$facet, "does not exist in the sample data. Faceting will be skipped."))
    }
}

# Save to output file using PDF device
print(paste("Saving plot to:", opt$output))
pdf(file = opt$output, width = 10, height = 8)
print(p)
dev.off()
