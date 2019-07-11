#!/usr/bin/env Rscript
#Normalises a SingleCellExperiment object

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

# Input from Loom format

scle <- import(opt$input_loom, format='loom', type='SingleCellLoomExperiment')
print(paste("Normalising...."))

#Normalise
scle <- normalize(scle, exprs_values = 1)

print(paste("Finished normalising"))

# Output to a Loom file
if (file.exists(opt$output_loom)) {
  file.remove(opt$output_loom)
}
export(scle, opt$output_loom, format='loom')
