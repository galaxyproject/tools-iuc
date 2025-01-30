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
    ),
    make_option(c("--nolines"),
        type = "logical", default = FALSE,
        help = "Remove borders (lines) around bars (TRUE/FALSE)"
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

# Debug: Check available taxonomic ranks
print("Available taxonomic ranks:")
print(colnames(tax_table(physeq)))

# Handle missing or unassigned taxa for all ranks
if (opt$keepNonAssigned) {
    # Replace NA or empty values with 'Not Assigned' for all ranks
    tax_df[is.na(tax_df) | tax_df == ""] <- "Not Assigned"
    tax_table(physeq) <- as.matrix(tax_df)
}

# Filter to top X taxa if requested
if (!is.null(opt$topX) && opt$topX != "") {
    topX <- as.numeric(opt$topX)
    if (is.na(topX) || topX <= 0) {
        stop("Error: topX should be a positive number.")
    }

    tax_rank <- opt$fill
    if (!tax_rank %in% colnames(tax_df)) {
        stop(paste("Error: Tax rank", tax_rank, "not found in tax_table."))
    }

    physeq_agg <- tax_glom(physeq, taxrank = tax_rank)
    taxa_abundance <- taxa_sums(physeq_agg)
    tax_table_agg <- tax_table(physeq_agg)
    taxa_abundance_by_rank <- tapply(taxa_abundance, tax_table_agg[, tax_rank], sum)
    top_taxa <- names(sort(taxa_abundance_by_rank, decreasing = TRUE))[1:topX]

    print("Top taxa:")
    print(top_taxa)

    otus_in_top_taxa <- rownames(tax_table_agg)[tax_table_agg[, tax_rank] %in% top_taxa]

    if (opt$keepOthers) {
        tax_table(physeq_agg)[, tax_rank][!rownames(tax_table_agg) %in% otus_in_top_taxa] <- "Others"
        physeq <- physeq_agg
    } else {
        physeq <- prune_taxa(otus_in_top_taxa, physeq_agg)
    }
}

# Generate bar plot
if (!is.null(opt$x) && opt$x != "") {
    p <- plot_bar(physeq, x = opt$x, fill = opt$fill) +
        geom_bar(aes(fill = !!sym(opt$fill)),
            stat = "identity", position = "stack",
            color = ifelse(opt$nolines, NA, "black")
        )
} else {
    p <- plot_bar(physeq, fill = opt$fill) +
        geom_bar(aes(fill = !!sym(opt$fill)),
            stat = "identity", position = "stack",
            color = ifelse(opt$nolines, NA, "black")
        )
}

# Optional: Add faceting if specified
if (!is.null(opt$facet) && opt$facet != "") {
    sample_vars <- colnames(sample_data(physeq))
    if (opt$facet %in% sample_vars) {
        p <- p + facet_wrap(as.formula(paste("~", opt$facet)))
    } else {
        warning(paste("Facet variable", opt$facet, "not found in sample data. Skipping faceting."))
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
