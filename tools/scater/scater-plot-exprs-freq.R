#!/usr/bin/env Rscript

#Plots mean expression vs % of expressing cells and provides information as to the number of genes expressed in 50% and 25% of cells.
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
    c("-o", "--output-plot-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path of the PDF output file to save plot to."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_loom', 'output_plot_file'))

# Check parameter values

if ( ! file.exists(opt$input_loom)){
  stop((paste('File', opt$input_loom, 'does not exist')))
}


# Input from Loom format

scle <- import(opt$input_loom, format='loom', type='SingleCellLoomExperiment')

#produce and save the scatter plot of reads vs genes
plot <- plotExprsFreqVsMean(scle, controls = "is_feature_control_MT")
ggsave(opt$output_plot_file, plot, device="pdf")
