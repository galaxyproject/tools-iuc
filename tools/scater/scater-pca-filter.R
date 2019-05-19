#!/usr/bin/env Rscript
#Filters a SingleCellExperiment object, using PCA on the following metrics:
# "pct_counts_top_100_features"
# "total_features"
# "pct_counts_feature_control"
# "total_features_feature_control"
# "log10_total_counts_endogenous"
# "log10_total_counts_feature_control"

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(mvoutlier))


# parse options
option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A serialized SingleCellExperiment object file in RDS format."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized SingleCellExperiment object."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

# Input from serialized R object

sce <- readRDS(opt$input_object_file)
print(paste("Starting with", ncol(sce), "cells and", nrow(sce), "features."))

#run PCA on data and detect outliers
sce <- runPCA(sce, use_coldata = TRUE, detect_outliers = TRUE, return_sce = TRUE)

#Identify and return non-outliers
sce$use <- !sce$outlier
sce <- sce[, colData(sce)$use]

print(paste("Ending with", ncol(sce), "cells and", nrow(sce), "features."))


# Output to a serialized R object
saveRDS(sce, file = opt$output_object_file)
