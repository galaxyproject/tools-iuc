# Load all required libraries
library(CEMiTool, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)
# setup R error handling to go to stderr
options(
  show.error.messages = FALSE,
  error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, FALSE)
  }
)

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

################################################################################
### Input Processing
################################################################################

# Collect arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "expressionMatrix", "M", 1, "character",
  "sampleAnnotation", "A", 2, "character",
  "pathwayList", "P", 2, "character",
  "interactions", "I", 2, "character",
  "filter","f", 1, "logical",
  "filter_pval", "i", 1, "numeric",
  "apply_vst", "a", 1, "logical",
  "n_genes", "n", 1, "integer",
  "eps", "e", 1, "numeric",
  "cor_method", "c", 1, "character",
  "cor_function", "y", 1, "character",
  "network_type", "x", 1, "character",
  "tom_type", "t", 1, "character",
  "merge_similar", "m", 1, "logical",
  "rank_method", "r", 1, "character",
  "min_ngen", "g", 1, "integer",
  "diss_thresh", "d", 1, "numeric",
  "center_func", "h", 1, "character",
  "ora_pval", "o", 1, "numeric",
  "gsea_scale", "l", 1, "logical",
  "gsea_min_size", "w", 1, "integer",
  "gsea_max_size", "z", 1, "integer"),
byrow = TRUE, ncol = 4
)

opt <- getopt(spec)
counts <- read.table(
  opt$expressionMatrix,
  header = TRUE,
  sep = "\t",
  strip.white = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Run CEMiTool

if (is.null(opt$sampleAnnotation)) {
  cem <- cemitool(
    counts,
    filter = opt$filter,
    filter_pval = opt$filter_pval,
    apply_vst = opt$apply_vst,
    n_genes = opt$n_genes,
    eps = opt$eps,
    cor_method = opt$cor_method,
    cor_function = opt$cor_function,
    network_type = opt$network_type,
    tom_type = opt$tom_type,
    merge_similar = opt$merge_similar,
    rank_method = opt$rank_method,
    min_ngen = opt$min_ngen,
    diss_thresh = opt$diss_thresh,
    center_func = opt$center_func,
    verbose = TRUE,
    ora_pval = opt$ora_pval
  )
} else {
  annotation <- read.table(
    opt$sampleAnnotation,
    header = TRUE,
    sep = "\t",
    strip.white = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cem <- cemitool(
    counts,
    annotation,
    filter = opt$filter,
    filter_pval = opt$filter_pval,
    apply_vst = opt$apply_vst,
    n_genes = opt$n_genes,
    eps = opt$eps,
    cor_method = opt$cor_method,
    cor_function = opt$cor_function,
    network_type = opt$network_type,
    tom_type = opt$tom_type,
    merge_similar = opt$merge_similar,
    rank_method = opt$rank_method,
    min_ngen = opt$min_ngen,
    diss_thresh = opt$diss_thresh,
    center_func = opt$center_func,
    verbose = TRUE,
    ora_pval = opt$ora_pval,
    gsea_scale = opt$gsea_scale,
    gsea_min_size = opt$gsea_min_size,
    gsea_max_size = opt$gsea_max_size,
    sample_name_column = "SampleName",
    class_column = "Class",
    order_by_class = TRUE
  )
  cem <- mod_gsea(cem)
  cem <- plot_gsea(cem)
}

if (!is.null(opt$pathwayList)) {
  ##print(opt$pathwayList)
  gmt_in <- read_gmt(opt$pathwayList)
  cem <- mod_ora(cem, gmt_in)
  cem <- plot_ora(cem)
}

if (!is.null(opt$interactions)) {
  interactions <- read.table(
    opt$interactions,
    header = TRUE,
    sep = "\t",
    strip.white = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ##print(opt$pathwayList)
  interactions_data(cem) <- interactions # add interactions
  cem <- plot_interactions(cem)
}

## Write analysis results into files
write_files(cem,
            directory = "./Tables",
            force = TRUE)

generate_report(cem)

save_plots(cem,
           value = "all",
           directory = "./Plots",
           force = TRUE)
