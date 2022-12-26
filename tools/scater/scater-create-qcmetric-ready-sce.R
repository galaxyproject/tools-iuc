#!/usr/bin/env Rscript
#Creates a SingleCellExperiment object, which scater's calculateQCMetrics already applied

library(optparse)
library(workflowscriptscommon)
library(scater)
library(LoomExperiment)

# parse options
#SCE-specific options
option_list <- list(
  make_option(
    c("-a", "--counts"),
    action = "store",
    default = NA,
    type = "character",
    help = "A tab-delimited expression matrix. The first column of all files is assumed to be feature names and the first row is assumed to be sample names."
  ),
  make_option(
    c("-r", "--row-data"),
    action = "store",
    default = NULL,
    type = "character",
    help = "Path to TSV (tab-delimited) format file describing the features. Row names from the expression matrix (-a), if present, become the row names of the SingleCellExperiment."
  ),
  make_option(
    c("-c", "--col-data"),
    action = "store",
    default = NULL,
    type = "character",
    help = "Path to TSV format file describing the samples (annotation). The number of rows (samples) must equal the number of columns in the expression matrix."
  ),
  #The scater-specific options
  make_option(
    c("-f", "--mt-controls"),
    action = "store",
    default = NULL,
    type = "character",
    help = "Path to file containing a list of the mitochondrial control genes"
  ),
  make_option(
    c("-p", "--ercc-controls"),
    action = "store",
    default = NULL,
    type = "character",
    help = "Path to file containing a list of the ERCC controls"
  ),
  make_option(
    c("-l", "--cell-controls"),
    action = "store",
    default = NULL,
    type = "character",
    help = "Path to file (one cell per line) to be used to derive a vector of cell (sample) names used to identify cell controls (for example, blank wells or bulk controls)."
  ),
  make_option(
    c("-o", "--output-loom"),
    action = "store",
    default = NA,
    type = "character",
    help = "File name in which to store the SingleCellExperiment object in Loom format."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c("counts", "output_loom"))

# Read the expression matrix

counts <- wsc_split_string(opt$counts)
reads <- read.table(counts)

# Read row and column annotations

rowdata <- opt$row_data

if (! is.null(opt$row_data)) {
  rowdata <- read.delim(opt$row_data)
}

coldata <- opt$col_data

if (! is.null(opt$col_data)) {
  coldata <- read.delim(opt$col_data)
}

# Now build the object

sce <- SingleCellLoomExperiment(assays = list(counts = as.matrix(reads)), colData = coldata)
#Scater options

# Check feature_controls (only mitochondrial and ERCC used for now)

if (! is.null(opt$mt_controls)) {
  if (! file.exists(opt$mt_controls)) {
    stop((paste("Supplied feature_controls file", opt$mt_controls, "does not exist")))
  } else {
    mts <- readLines(opt$mt_controls)
  }
} else {
  mts <- NULL
}

if (! is.null(opt$ercc_controls)) {
  if (! file.exists(opt$ercc_controls)) {
    stop((paste("Supplied feature_controls file", opt$ercc_controls, "does not exist")))
  } else {
    ercc_controls <- readLines(opt$ercc_controls)
  }
} else {
  ercc_controls <- NULL
}

# Check cell_controls

if (! is.null(opt$cell_controls)) {
  if (! file.exists(opt$cell_controls)) {
    stop((paste("Supplied feature_controls file", opt$cell_controls, "does not exist")))
  } else {
    cell_controls <- readLines(opt$cell_controls)
  }
} else {
  cell_controls <- NULL
}

# calculate QCMs

sce <- addPerCellQC(sce, subsets = list(Mito = mts, ERCC = ercc_controls, cell_controls = cell_controls))

# Output to a Loom file
if (file.exists(opt$output_loom)) {
  file.remove(opt$output_loom)
}
export(sce, opt$output_loom, format = "loom")
