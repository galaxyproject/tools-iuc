#!/usr/bin/env Rscript
# Manually filter SingleCellExperiment with user-defined parameters

# Load optparse we need to check inputs
library(optparse)
library(workflowscriptscommon)
library(LoomExperiment)
library(scater)

# parse options
option_list = list(
  make_option(
    c("-i", "--input-loom"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A SingleCellExperiment object file in Loom format."
  ),
  make_option(
    c("-d", "--detection-limit"),
    action = "store",
    default = 0,
    type = 'numeric',
    help = "Numeric scalar providing the value above which observations are deemed to be expressed"
  ),
  make_option(
    c("-l", "--library-size"),
    action = "store",
    default = 0,
    type = 'numeric',
    help = "Minimum library size (mapped reads) to filter cells on"
  ),
  make_option(
    c("-e", "--expressed-genes"),
    action = "store",
    default = 0,
    type = 'numeric',
    help = "Minimum number of expressed genes to filter cells on"
  ),
  make_option(
    c("-m", "--percent-counts-MT"),
    action = "store",
    default = 100,
    type = 'numeric',
    help = "Maximum % of mitochondrial genes expressed per cell. Cells that exceed this value will be filtered out."
  ),
  make_option(
    c("-o", "--output-loom"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store the SingleCellExperiment object in Loom format."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_loom', 'output_loom'))

# Check parameter values

if ( ! file.exists(opt$input_loom)){
  stop((paste('File', opt$input_loom, 'does not exist')))
}

# Filter out unexpressed features

scle <- import(opt$input_loom, format='loom', type='SingleCellLoomExperiment')
print(paste("Starting with", ncol(scle), "cells and", nrow(scle), "features."))

# Create a logical vector of features that are expressed (above detection_limit)
feature_expressed <- nexprs(scle, detection_limit = opt$detection_limit, exprs_values = 1, byrow=TRUE) > 0
scle <- scle[feature_expressed, ]

print(paste("After filtering out unexpressed features: ", ncol(scle), "cells and", nrow(scle), "features."))

# Filter low library sizes
to_keep <- scle$total_counts > opt$library_size
scle <- scle[, to_keep]

print(paste("After filtering out low library counts: ", ncol(scle), "cells and", nrow(scle), "features."))


# Filter low expressed genes
to_keep <- scle$total_features > opt$expressed_genes
scle <- scle[, to_keep]

print(paste("After filtering out low expressed: ", ncol(scle), "cells and", nrow(scle), "features."))


# Filter out high MT counts
to_keep <- scle$pct_counts_MT < opt$percent_counts_MT
scle <- scle[, to_keep]

print(paste("After filtering out high MT gene counts: ", ncol(scle), "cells and", nrow(scle), "features."))

# Output to a Loom file
if (file.exists(opt$output_loom)) {
  file.remove(opt$output_loom)
}
export(scle, opt$output_loom, format='loom')
