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

parser$add_argument("--factorLevel", action = "append", required = TRUE)
parser$add_argument("--factorLevel_counts",
                    action = "append",
                    required = TRUE)
parser$add_argument("--factorLevel_n", action = "append",  required = TRUE)
parser$add_argument("--cores",  type = "integer", required = TRUE)
parser$add_argument("--normalize", action = "store_true", required = FALSE)
parser$add_argument("--nbins", type = "integer", required = TRUE)
parser$add_argument("--lwr", type = "numeric", required = TRUE)
parser$add_argument("--upr", type = "numeric", required = TRUE)

args <- parser$parse_args()

all_files <- args$factorLevel_counts

conditions <- c()
for (x in 1:length(args$factorLevel)) {
  temp <- append(conditions, rep(args$factorLevel[[x]]))
  conditions <- temp
}

sample_names <- all_files %>%
  str_replace(pattern = "\\.tab", "")

design <-
  data.frame(list(
    sample = sample_names,
    condition = conditions,
    path = all_files
  ))
so <- sleuth_prep(design,
                  cores = args$cores,
                  normalize = args$normalize)

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
