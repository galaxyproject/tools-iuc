#!/usr/bin/env Rscript

## Redirect R error handling to stderr.
options(show.error.messages = FALSE,
        error = function() {
          cat(geterrmessage(), file = stderr())
          q("no", 1, FALSE)
        })

## Avoid crashing Galaxy with a UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

args <- commandArgs(TRUE)
if (length(args) == 0) {
  stop("Arguments missing for Rscrpit", call. = FALSE)
} else {
  # percentage of identity
  id_threshold <- as.numeric(args[3])
  # get input data (matrix)
  data <- read.csv(args[1], header = FALSE, sep = ",", row.names = 1)
  # remove last 2 columns
  data_length <- length(data)
  # create matrix
  mat <- as.matrix(data[, 1:data_length], fill = TRUE)
  # create coordinate matrix
  d <- as.dist(1 - mat)
  # create tree
  hc <- hclust(d, method = "single")
  # assign otu based on identity value
  otu <- cutree(hc, h = -id_threshold)
  # group contigs by otu
  # Print results to output file
  output <- args[2]
  # unique is used to know the number of different otu
  for (i in unique(otu)){
    # retrieve contigs belonging to the same otu
    clust <- which(otu == i)
    # write otu number and number of contigs in this otu
    cat(
        paste("OTU_", i, ",", length(clust), ",", sep = ""),
        file = output, append = TRUE)
    for (n in names(clust)){
      # write contigs name
      cat(paste(gsub(" ", "", n), ",", sep = ""), file = output, append = TRUE)
    }
    cat("\n", sep = "", file = output, append = TRUE)
  }
}
