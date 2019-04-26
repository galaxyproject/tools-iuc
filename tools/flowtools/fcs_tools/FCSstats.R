#!/usr/bin/Rscript --vanilla
# FCS Summary Statistic Module for Galaxy
# FlowCore
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#

library(flowCore)

getMarkerNames <- function(input, output) {
  fcs <- read.FCS(input, transformation=F)

  ## marker names
  channels <- colnames(fcs)
  markers <- as.vector(pData(parameters(fcs))$desc)
  df <- data.frame(channels, markers)
  fcs_summary <- capture.output(summary(fcs))
  fcs_dim <- capture.output(dim(fcs))
  fcs_markers <- capture.output(df)


  sink(output)
  cat(fcs_dim, sep="\n")
  cat("\n\n=========================\n")
  cat("==     FCS SUMMARY     ==\n")
  cat("=========================\n")
  cat(fcs_summary, sep="\n")
  cat("\n\n=========================\n")
  cat("==    MARKERS IN FCS   ==\n")
  cat("=========================\n")
  cat(fcs_markers, sep="\n")
  sink()
}

checkFCS <- function(input_file, output_file) {
  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid <- isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })

  if (isValid) {
    getMarkerNames(input_file, output_file)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
checkFCS(args[1], args[2])
