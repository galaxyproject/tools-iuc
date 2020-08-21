#!/usr/bin/Rscript --vanilla
# Cell Ontology Module for Galaxy
# FlowCL
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#

library(flowCL)
library(base)

getOntology <- function(output_file, markers) {
  res <- flowCL(markers, ResetArch = TRUE)
  if (length(res) == 6) {
    report <- capture.output(res$Table)
    sink(output_file)
    cat(report, sep = "\n")
    sink()
  }
}

args <- commandArgs(trailingOnly = TRUE)
markers <- paste(args[2:length(args)], collapse="")
getOntology(args[1], markers)
