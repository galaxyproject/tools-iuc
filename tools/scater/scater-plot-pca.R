#!/usr/bin/env Rscript

# Creates a PCA plot of a normalised SingleCellExperiment object.

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
    c("-c", "--colour-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to colour PCA plot points by. The values represented in this options should be categorical"
  ),
  make_option(
    c("-s", "--size-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to size PCA plot points by. The values represented in this options should be numerical and not categorical"
  ),
  make_option(
    c("-p", "--shape-by"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Feature (from annotation file) to shape PCA plot points by. The values represented in this options should be categorical"
  ),
  make_option(
    c("-o", "--output-plot-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path of the PDF output file to save plot to."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_plot_file'))
# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}


# Input from serialized R object

sce <- readRDS(opt$input_object_file)
sce <- normalize(sce)
sce <- runPCA(sce)
plot <- plotReducedDim(sce, "PCA", colour_by = opt$colour_by, size_by = opt$size_by, shape_by = opt$shape_by)
#do the scatter plot of reads vs genes

ggsave(opt$output_plot_file, plot, device="pdf")
