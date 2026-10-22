#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

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
        help = "Variable for fill color (e.g., 'Genus', 'Order'). Use 'ASV' as argument to show each OTU/ASV."
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
    make_option(c("--normalize_x"),
        action = "store_true", dest = "normalize_x", default = FALSE,
        help = "Normalize x groups to sum up to 100%"
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
        help = "Output format (e.g., 'pdf', 'png', 'jpg')"
    ),
    make_option(c("--nolines"),
        action = "store_true", dest = "nolines", default = FALSE,
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

if (is.null(opt$fill) || opt$fill == "") {
    print(paste("No fill chosen using ASV"))
    opt$fill <- "ASV"
}

# Load phyloseq object
print(paste("Trying to read:", opt$input))
physeq <- readRDS(opt$input)

## Allow to use OTU as tax group
# Extract rownames (taxids) from the tax_table and add them as a new column
taxids <- rownames(tax_table(physeq))

# Get the number of columns in the tax_table
num_columns <- ncol(tax_table(physeq))

# Add the taxids as a new last column in the tax_table
tax_table(physeq) <- cbind(tax_table(physeq), taxid = taxids)

# Rename the last column to 'ASV' / OTU does conflict with phyloseq logic
colnames(tax_table(physeq))[num_columns + 1] <- "ASV"

# Normalize to relative abundances if requested
if (opt$normalize) {
    print("Normalizing abundances to sum to 100%...")
    physeq <- transform_sample_counts(physeq, function(x) 100 * x / sum(x))
}

# Debug: Check available taxonomic ranks

tax_ranks <- colnames(tax_table(physeq))
sample_vars <- colnames(sample_data(physeq))

print("Available taxonomic ranks:")
print(tax_ranks)

print("Available metadata:")
print(sample_vars)

# Handle missing or unassigned taxa for all ranks
if (opt$keepNonAssigned) {
    # Replace NA or empty values with 'Not Assigned' for all ranks

    for (rank in tax_ranks) {
        if (rank %in% tax_ranks) {
            # replace NA or empty values with 'Not Assigned'
            tax_table(physeq)[, rank][is.na(tax_table(physeq)[, rank])] <- "Not Assigned"
        }
    }
}

# Filter to top X taxa if requested
if (!is.null(opt$topX) && opt$topX != "") {
    topX <- as.numeric(opt$topX)
    if (is.na(topX) || topX <= 0) {
        stop("Error: topX should be a positive number.")
    }

    tax_rank <- opt$fill
    if (!tax_rank %in% colnames(tax_table(physeq))) {
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

    # Group non-top OTUs as 'Others' if requested
    if (opt$keepOthers) {
        # Update the tax_table to assign 'Others' to non-top taxa
        tax_table(physeq_agg)[, tax_rank][!rownames(tax_table_agg) %in% otus_in_top_taxa] <- "Others"
        physeq <- physeq_agg
    } else {
        physeq <- prune_taxa(otus_in_top_taxa, physeq_agg)
    }
}


# normalize x groups if needed
if (opt$x %in% sample_vars) {
    if (opt$normalize_x && !is.null(opt$x) && opt$x != "") {
        physeq_agg <- merge_samples(physeq, opt$x)

        physeq <- transform_sample_counts(physeq_agg, function(x) (x / sum(x) * 100))
        opt$x <- NULL # set to Null since we do not need x for downstream now
        opt$facet <- NULL # set to Null since facetting does not work with normalize x
        warning(paste("normalize x does not work with facetting"))
    }
} else {
    warning(paste("x", opt$x, "not found in sample data. Skipping normalize_x."))
}


# Check if the facet variable is valid and exists
facet_var <- NULL
if (!is.null(opt$facet) && opt$facet != "") {
    if (opt$facet %in% sample_vars || opt$facet %in% tax_ranks) {
        facet_var <- opt$facet # Store facet variable for later
    } else {
        warning(paste("Facet variable", opt$facet, "not found in sample data or tax ranks. Skipping faceting."))
    }
}

# Determine if faceting is needed
facet_formula <- if (!is.null(facet_var)) as.formula(paste("~", facet_var)) else NULL

# Define color based on the `nolines` option
plot_color <- ifelse(opt$nolines, NA, "black")

# Generate bar plot
if (!is.null(opt$x) && opt$x != "") {
    p <- plot_bar(physeq,
        x = opt$x,
        fill = opt$fill
    ) + facet_wrap(facet_formula, scales = "free_x") +
        geom_bar(
            stat = "identity",
            position = "stack",
            aes(fill = !!sym(opt$fill)),
            color = plot_color
        )
} else {
    p <- plot_bar(physeq,
        fill = opt$fill
    ) + facet_wrap(facet_formula, scales = "free_x") +
        geom_bar(
            stat = "identity",
            position = "stack",
            aes(fill = !!sym(opt$fill)),
            color = plot_color
        )
}


# Reorder fill levels to ensure "Not Assigned" and "Others" are at the bottom if they exist
fill_values <- unique(p$data[[opt$fill]]) # Get unique fill values
new_levels <- setdiff(fill_values, c("Not Assigned", "Others")) # Exclude "Not Assigned" and "Others"

if ("Not Assigned" %in% fill_values) {
    new_levels <- c("Not Assigned", new_levels) # Place "Not Assigned" at the bottom if it exists
}
if ("Others" %in% fill_values) {
    new_levels <- c("Others", new_levels) # Place "Others" at the bottom if it exists
}

# Apply the new levels to the fill variable in the plot data
p$data[[opt$fill]] <- factor(p$data[[opt$fill]], levels = new_levels)

# Save to output file
ggsave(
    filename = opt$output,
    plot = p,
    width = opt$width,
    height = opt$height,
    device = opt$device
)
