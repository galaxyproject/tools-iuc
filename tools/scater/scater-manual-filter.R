#!/usr/bin/env Rscript
# Manually filter SingleCellExperiment with user-defined parameters

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(scater))

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
    c("-m", "--percent-counts-MT"),
    action = "store",
    default = 100,
    type = 'numeric',
    help = "Maximum % of mitochondrial genes expressed per cell. Cells that exceed this value will be filtered out."
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

# Filter out unexpressed features

sce <- readRDS(opt$input_object_file)
print(paste("Starting with", ncol(sce), "cells and", nrow(sce), "features."))

#create a logical vector of features that are expressed (obove detection_limit)
feature_expressed <- nexprs(sce, detection_limit = opt$detection_limit, byrow=TRUE) > 0
sce <- sce[feature_expressed, ]

print(paste("After filtering out unexpressed features: ", ncol(sce), "cells and", nrow(sce), "features."))

#Filter low library sizes
sce$use <- sce$total_counts > opt$library_size
sce <- sce[, colData(sce)$use]

print(paste("After filtering out low library counts: ", ncol(sce), "cells and", nrow(sce), "features."))

#Filter out high MT counts

sce$use <- sce$pct_counts_MT < opt$percent_counts_MT
sce <- sce[, colData(sce)$use]

print(paste("After filtering out high MT gene counts: ", ncol(sce), "cells and", nrow(sce), "features."))


# Output to a serialized R object
saveRDS(sce, file = opt$output_object_file)
