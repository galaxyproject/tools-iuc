#!/usr/bin/env Rscript

#Plots mean expression vs % of expressing cells and provides information as to the number of genes expressed in 50% and 25% of cells.
# Load optparse we need to check inputs

suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(optparse))


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
    c("-o", "--output-plot-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A jpg file to save plot to, e.g Reads.plot."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_plot_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}


# Input from serialized R object

sce <- readRDS(opt$input_object_file)

#produce and save the scatter plot of reads vs genes
plot <- plotExprsFreqVsMean(sce,controls = "is_feature_control_MT")
ggsave(opt$output_plot_file, plot, device="pdf")
