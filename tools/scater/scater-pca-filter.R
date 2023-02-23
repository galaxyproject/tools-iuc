#!/usr/bin/env Rscript
#Filters a SingleCellExperiment object, using PCA on the following metrics:
# "pct_counts_top_100_features"
# "total_features"
# "pct_counts_feature_control"
# "total_features_feature_control"
# "log10_total_counts_endogenous"
# "log10_total_counts_feature_control"

# Load optparse we need to check inputs
library(optparse)
library(workflowscriptscommon)
library(LoomExperiment)
library(scater)
library(robustbase)

# parse options
option_list <- list(
  make_option(
    c("-i", "--input-loom"),
    action = "store",
    default = NA,
    type = "character",
    help = "A SingleCellExperiment object file in Loom format."
  ),
  make_option(
    c("-o", "--output-loom"),
    action = "store",
    default = NA,
    type = "character",
    help = "File name in which to store the SingleCellExperiment object in Loom format."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c("input_loom", "output_loom"))

# Check parameter values

if (! file.exists(opt$input_loom)) {
  stop((paste("File", opt$input_loom, "does not exist")))
}

# Filter out unexpressed features

sce <- import(opt$input_loom, format = "loom", type = "SingleCellLoomExperiment")

print(paste("Starting with", ncol(sce), "cells"))

sce <- runColDataPCA(sce, outliers = TRUE, variables = list("sum", "detected", "subsets_Mito_percent"))
sce$use <- !sce$outlier
sce <- sce[, colData(sce)$use]
print(paste("Ending with", ncol(sce), "cells"))


# Output to a Loom file
if (file.exists(opt$output_loom)) {
  file.remove(opt$output_loom)
}
export(sce, opt$output_loom, format = "loom")
