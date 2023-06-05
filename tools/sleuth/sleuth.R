library(sleuth,
        quietly = TRUE,
        warn.conflicts = FALSE)
library(annotables, quietly = TRUE, warn.conflicts = FALSE)
library(argparse, quietly = TRUE, warn.conflicts = FALSE)
library(tidyverse)


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
parser <- ArgumentParser(description = "Sleuth R script")

parser$add_argument("--factorLevel", action = "append", required = FALSE)
parser$add_argument("--factorLevel_counts",
                    action = "append",
                    required = FALSE)
parser$add_argument("--factorLevel_n", action = "append",  required = FALSE)
parser$add_argument("--cores",  type = "integer", required = FALSE)
parser$add_argument("--normalize", action = "store_true", required = FALSE)
parser$add_argument("--nbins", type = "integer", required = FALSE)
parser$add_argument("--lwr", type = "numeric", required = FALSE)
parser$add_argument("--upr", type = "numeric", required = FALSE)
parser$add_argument("--metadata_file",
                    action = "append",
                    required = FALSE)
parser$add_argument("--experiment_design", required = FALSE)

args <- parser$parse_args()

if (args$experiment_design == "complex") {
  ## Complex experiment design
  ############################
  
  s2c  <-
    read.table(file = args$metadata_file,
               header = TRUE,
               sep = "\t")
  paths <- c()
  for (x in s2c$data_filename) {
    paths <- c(paths, paste('./kallisto_outputs/', x, sep = ""))
  }
  for (f in paths) {
    file.rename(f, gsub(".fastq.*", ".h5", f))
  }
  s2c$path <- gsub(".fastq.*", ".h5", paths)
  
  so <- sleuth_prep(s2c, full_model = ~ condition, num_cores = 1)
  so <- sleuth_fit(so)
  
} else {
  ## Simple experiment design
  ###########################
  
  conditions <- c()
  for (x in seq_along(args$factorLevel)) {
    temp <- append(conditions, rep(args$factorLevel[[x]]))
    conditions <- temp
  }
  
  sample_names <-
    gsub(".fastq.+", "", basename(args$factorLevel_counts))
  
  design <-
    data.frame(list(
      sample = sample_names,
      condition = conditions,
      path = args$factorLevel_counts
    ))
  so <- sleuth_prep(design,
                    cores = args$cores,
                    normalize = args$normalize)
}

so <- sleuth_fit(
  so,
  ~ condition,
  "full",
  n_bins = args$nbins,
  lwr = args$lwr,
  upr = args$upr
)

so <- sleuth_fit(
  so,
  ~ 1,
  "reduced",
  n_bins = args$nbins,
  lwr = args$lwr,
  upr = args$upr
)

so <- sleuth_lrt(so, "reduced", "full")

sleuth_table <-
  sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)

write.table(
  sleuth_table,
  file = "sleuth_table.tab",
  quote = FALSE,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)


outputFile <- file.path(getwd(), "pca_plot.pdf")
pdf(file = outputFile,
    height = 6,
    width = 9)
plot_pca(so, color_by = "condition")
dev.off()

outputFile <- file.path(getwd(), "group_density.pdf")
pdf(file = outputFile,
    height = 6,
    width = 9)
plot_group_density(
  so,
  use_filtered = TRUE,
  units = "est_counts",
  trans = "log",
  grouping = setdiff(colnames(so$sample_to_covariates),
                     "sample"),
  offset = 1
)
dev.off()