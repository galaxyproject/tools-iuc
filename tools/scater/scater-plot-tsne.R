#!/usr/bin/env Rscript

# Creates a t-SNE plot of a normalised SingleCellExperiment object.

# Load optparse we need to check inputs

library(optparse)
library(workflowscriptscommon)
library(LoomExperiment)
library(scater)
library(Rtsne)

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
    c("-c", "--colour-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to colour t-SNE plot points by. The values represented in this options should be categorical"
  ),
  make_option(
    c("-s", "--size-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to size t-SNE plot points by. The values represented in this options should be numerical and not categorical"
  ),
  make_option(
    c("-p", "--shape-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to shape t-SNE plot points by. The values represented in this options should be categorical"
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
scle <- normalize(scle, exprs_values = 1)
scle <- runTSNE(scle, perplexity=10)
plot <- plotTSNE(scle, colour_by = opt$colour_by, size_by = opt$size_by, shape_by = opt$shape_by)


ggsave(opt$output_plot_file, plot, device="pdf")
