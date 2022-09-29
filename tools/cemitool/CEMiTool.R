# Load all required libraries
library(CEMiTool, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(getopt, quietly = TRUE, warn.conflicts = FALSE)
# setup R error handling to go to stderr
options(show.error.messages = F, error = function() {
  cat(geterrmessage(), file = stderr())
  q("no", 1, F)
})

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
  "interactions", "I", 2, "character"),
byrow = TRUE, ncol = 4
)

opt <- getopt(spec)
counts <- read.table(opt$expressionMatrix,
                      header = TRUE, 
                      sep = "\t", 
                      strip.white = TRUE,
                      stringsAsFactors = FALSE, 
                      check.names = FALSE)


# Run CEMiTool

if (is.null(opt$sampleAnnotation)) {
  cem <- cemitool(counts,
                  filter = TRUE,
                  filter_pval = 0.1,
                  apply_vst = FALSE,
                  n_genes=1000,
                  eps = 0.1,
                  cor_method = c("pearson", "spearman"),
                  cor_function = "cor",
                  network_type = "unsigned",
                  tom_type = "signed",
                  set_beta = NULL,
                  force_beta = FALSE,
                  merge_similar = TRUE,
                  rank_method = "mean",
                  min_ngen = 30,
                  diss_thresh = 0.8,
                  center_func = "mean",
                  directed = FALSE,
                  verbose = TRUE,
                  ora_pval = 0.05)
} else {
  annotation <- read.table(opt$sampleAnnotation,
                           header = TRUE, 
                           sep = "\t", 
                           strip.white = TRUE, 
                           stringsAsFactors = FALSE, 
                           check.names = FALSE)
  cem <- cemitool(counts,
                  annotation,
                  gsea_scale = TRUE,
                  gsea_min_size = 15,
                  gsea_max_size = 1000,
                  sample_name_column = "SampleName",
                  class_column = "Class",
                  order_by_class = TRUE,
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
  interactions <- read.table(opt$interactions,
                           header = TRUE, 
                           sep = "\t", 
                           strip.white = TRUE, 
                           stringsAsFactors = FALSE, 
                           check.names = FALSE)
  ##print(opt$pathwayList)
  interactions_data(cem) <- interactions # add interactions
  cem <- plot_interactions(cem)
}

## Write analysis results into files
write_files(cem, 
            directory="./Tables",
            force=TRUE)

generate_report(cem)

save_plots(cem,
           value="all",
           directory="./Plots",
           force=TRUE)

## Generate GSEA plot
##cem <- mod_gsea(cem)
##cem <- plot_gsea(cem)
## Save gsea plots


##data(expr0)
# 
# write.table(int_df, file = "interactions.tab",
#             sep = "\t", row.names = F)
# 
# n_genes <- 100
# # run cemitool
# cemitool(
#   expr0,
#   sample_annot,
#   gmt_in,
#   int_df,
#   filter = TRUE,
#   filter_pval = 0.1,
#   apply_vst = FALSE,
#   n_genes,
#   eps = 0.1,
#   cor_method = c("pearson", "spearman"),
#   cor_function = "cor",
#   network_type = "unsigned",
#   tom_type = "signed",
#   set_beta = NULL,
#   force_beta = FALSE,
#   sample_name_column = "SampleName",
#   class_column = "Class",
#   merge_similar = TRUE,
#   rank_method = "mean",
#   ora_pval = 0.05,
#   gsea_scale = TRUE,
#   gsea_min_size = 15,
#   gsea_max_size = 1000,
#   min_ngen = 30,
#   diss_thresh = 0.8,
#   plot = TRUE,
#   plot_diagnostics = TRUE,
#   order_by_class = TRUE,
#   center_func = "mean",
#   directed = FALSE,
#   verbose = FALSE
# )
# 
# cemitool(
#   expr0,
#   ##sample_annot,
#   filter = TRUE,
#   filter_pval = 0.1,
#   apply_vst = FALSE,
#   n_genes = 100,
#   eps = 0.1,
#   cor_method = "pearson",
#   cor_function = "cor",
#   network_type = "unsigned",
#   tom_type = "signed",
#   set_beta = NULL,
#   force_beta = FALSE,
#   sample_name_column = "SampleName",
#   class_column = "Class",
#   merge_similar = TRUE,
#   rank_method = "mean",
#   ora_pval = 0.05,
#   gsea_scale = TRUE,
#   gsea_min_size = 15,
#   gsea_max_size = 1000,
#   min_ngen = 30,
#   diss_thresh = 0.8,
#   plot = TRUE,
#   plot_diagnostics = TRUE,
#   order_by_class = TRUE,
#   center_func = "mean",
#   directed = FALSE,
#   verbose = FALSE
# )
# 
# ## Write analysis results into files
# write_files(cem, 
#             directory="./Tables",
#             force=TRUE)
# 
# ## Generate GSEA plot
# cem <- mod_gsea(cem)
# cem <- plot_gsea(cem)
# 
# 
# 
# 
# ## Save gsea plots
# save_plots(cem,
#            value="gsea",
#            directory="./Plots",
#            force=TRUE)
# 
# 
# 
# ## Save all plots
# save_plots(cem,
#            value="all",
#            directory="./Plots")
# 
# # create report as html document
# generate_report(cem, 
#                 directory="./Report",
#                 force=TRUE)
# 
# 
# 
# # save all plots
# save_plots(cem, "all", directory="./Plots")
# ## ------------------------------------------------------------------------
# ## 
# ## ------------------------------------------------------------------------