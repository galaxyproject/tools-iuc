# FCS Headers Module for Galaxy
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

getFCSMarkerNames <- function(input, output) {
  fcs <- read.FCS(input, transformation=F)
  ## marker names
  channels <- colnames(fcs)
  markers <- as.vector(pData(parameters(fcs))$desc)
  df <- data.frame(channels, markers)
  fcs_markers <- capture.output(df)

  write.table(df, output, sep="\t")
}

checkFCS <- function(input_file, output_file) {
  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid = isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })

  if (isValid) {
    getFCSMarkerNames(input_file, output_file)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
checkFCS(args[2], args[3])
