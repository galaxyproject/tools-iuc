#!/usr/bin/env Rscript
# Manually filter SingleCellExperiment with user-defined parameters

# Load optparse we need to check inputs
library(optparse)
library(workflowscriptscommon)
library(LoomExperiment)
library(scater)

# parse options
option_list <- list(
    make_option(
        c("-i", "--input-loom"),
        action = "store",
        default = NA,
        type = "character",
        help = "A SingleCellExperiment object file in Loom format"
    ),
    make_option(
        c("-l", "--library-size"),
        action = "store",
        default = 0,
        type = "numeric",
        help = "Minimum library size (mapped reads) to filter cells on"
    ),
    make_option(
        c("-m", "--percent-counts-MT"),
        action = "store",
        default = 100,
        type = "numeric",
        help = "Maximum % of mitochondrial genes expressed per cell. Cells that exceed this value will be filtered out"
    ),
    make_option(
        c("-f", "--expressed-features"),
        action = "store",
        default = 100,
        type = "numeric",
        help = "Remove cells that have less than the given number of expressed features"
    ),
    make_option(
        c("-d", "--detection-limit"),
        action = "store",
        default = 0,
        type = "numeric",
        help = "Number of reads mapped to a feature above which it to be deemed as expressed"
    ),
    make_option(
        c("-e", "--min-cells-expressed"),
        action = "store",
        default = 0,
        type = "numeric",
        help = "Remove features that occur in less than the given number of cells"
    ),
    make_option(
        c("-o", "--output-loom"),
        action = "store",
        default = NA,
        type = "character",
        help = "File name in which to store the SingleCellExperiment object in Loom format"
    )
)

opt <- wsc_parse_args(option_list, mandatory = c("input_loom", "output_loom"))

# Check parameter values

if (!file.exists(opt$input_loom)) {
    stop((paste("File", opt$input_loom, "does not exist")))
}

# Filter out unexpressed features

sce <- import(opt$input_loom, format = "loom", type = "SingleCellLoomExperiment")
print(paste("Starting with", ncol(sce), "cells and", nrow(sce), "features."))

# Filter out low quality cells

# Filter low library sizes
passing_total <- sce$total > opt$library_size
sce <- sce[, passing_total]
print(paste("After filtering out low library counts: ", ncol(sce), "cells and", nrow(sce), "features."))

# Filter out high MT counts
passing_mt_counts <- sce$subsets_Mito_percent < opt$percent_counts_MT
sce <- sce[, passing_mt_counts]
print(paste("After filtering out high MT gene counts: ", ncol(sce), "cells and", nrow(sce), "features."))

expr_features <- sce$detected > opt$expressed_features
sce <- sce[, expr_features]
print(paste("After filtering out cells with low feature counts: ", ncol(sce), "cells and", nrow(sce), "features."))

# Create a logical vector of features that are expressed (above detection_limit)
feature_expressed <- nexprs(sce, detection_limit = opt$detection_limit, byrow = TRUE) > opt$min_cells_expressed
sce <- sce[feature_expressed, ]
print(paste("After filtering out rare features: ", ncol(sce), "cells and", nrow(sce), "features."))

# Output to a Loom file
if (file.exists(opt$output_loom)) {
    file.remove(opt$output_loom)
}
export(sce, opt$output_loom, format = "loom")
