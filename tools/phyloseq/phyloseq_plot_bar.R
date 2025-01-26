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
    make_option(c("--keepNonAssigned"),
        action = "store_true", dest = "keepNonAssigned", default = FALSE,
        help = "Keep taxa labeled as 'Not Assigned' (optional)"
    ),
    make_option(c("--normalize"),
        action = "store_true", dest = "normalize", default = FALSE,
        help = "Normalize abundances to sum to 100% (optional)"
    ),
    make_option(c("--width"),
        action = "store", dest = "width", default = 10,
        type = "numeric", help = "Width of the output plot in inches"
    ),
    make_option(c("--height"),
        action = "store", dest = "height", default = 8,
        type = "numeric", help = "Height of the output plot in inches"
    ),
    make_option(c("--device"),
        action = "store", dest = "device", default = "pdf",
        help = "Output format (e.g., 'pdf', 'png', 'jpeg')"
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

# Normalize to relative abundances if requested
if (opt$normalize) {
    print("Normalizing abundances to sum to 100%...")
    physeq <- transform_sample_counts(physeq, function(x) 100 * x / sum(x))
}

if (opt$keepNonAssigned) {
    # Replace missing/NA taxa with "Not Assigned"
    tax_table(physeq) <- apply(tax_table(physeq), c(1, 2), function(x) ifelse(is.na(x) | x == "", "Not Assigned", x))
} else {
    # Filter out "Not Assigned" taxa
    physeq <- subset_taxa(physeq, !is.na(TaxaLevel) & TaxaLevel != "Not Assigned")
}
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

    # Summarize the abundance at each taxonomic rank (grouping by taxonomic name)
    tax_table_agg <- tax_table(physeq_agg)
    taxa_abundance_by_rank <- tapply(taxa_abundance, tax_table_agg[, tax_rank], sum)

    # Identify the top X taxa by summed abundance
    top_taxa <- names(sort(taxa_abundance_by_rank, decreasing = TRUE))[1:topX]

    print("Only plotting taxa in TopX taxa group:")
    print(top_taxa)

    # Get the OTUs corresponding to the top taxa
    otus_in_top_taxa <- rownames(tax_table_agg)[tax_table_agg[, tax_rank] %in% top_taxa]

    if (opt$keepOthers) {
        # Label taxa not in top_taxa as "Others"
        tax_table(physeq_agg)[, tax_rank][!rownames(tax_table(physeq_agg)) %in% otus_in_top_taxa] <- "Others"
        physeq <- physeq_agg
    } else {
        # Subset the phyloseq object to keep only the top X taxa
        physeq_filtered <- prune_taxa(otus_in_top_taxa, physeq_agg)
        physeq <- physeq_filtered
    }
}

# Generate bar plot
if (!is.null(opt$x) && opt$x != "") {
    p <- plot_bar(physeq, x = opt$x, fill = opt$fill) +
        geom_bar(aes(color = NULL, fill = !!sym(opt$fill)), stat = "identity", position = "stack", color = NA)
} else {
    p <- plot_bar(physeq, fill = opt$fill) +
        geom_bar(aes(color = NULL, fill = !!sym(opt$fill)), stat = "identity", position = "stack", color = NA)
}

# Only facet if the facet variable is provided and exists in the sample data
if (!is.null(opt$facet) && opt$facet != "") {
    if (opt$facet %in% sample_vars) {
        p <- p + facet_wrap(as.formula(paste("~", opt$facet)))
    } else {
        warning(paste("Facet variable", opt$facet, "does not exist in the sample data. Faceting will be skipped."))
    }
}


# Save to output file
ggsave(
    filename = opt$output,
    plot = p,
    width = opt$width,
    height = opt$height,
    device = opt$device
)
