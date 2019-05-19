#!/usr/bin/env Rscript
#Normalises a SingleCellExperiment object

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
print(paste("Normalising...."))

#Normalise
sce <- normalize(sce)

print(paste("Finished normalising"))

# Output to a serialized R object
saveRDS(sce, file = opt$output_object_file)
