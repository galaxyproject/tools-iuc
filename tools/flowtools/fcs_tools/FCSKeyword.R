#!/usr/bin/Rscript --vanilla
# ImmPort FCSKeywords
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Converts the FCS file to text without transformaton
# To run in R
# 1) library(flowCore)
# 2) source("FCSKeyword.R")
# 3) transformFCS("filename")
#
# Version 1.4.1
# March 2016 -- added lines to run directly from command line
#

library(flowCore)

#
# Starting function for processing a FCS file
#
extractKeywords <- function(input_file, keyword_file="", debug=FALSE) {
  #
  # Generate the file names for the output_file and keyword_file
  #
  pieces <- unlist(strsplit(input_file, .Platform$file.sep))
  filename <- pieces[length(pieces)]

  if (keyword_file == "") {
    filepieces <- unlist(strsplit(filename, '\\.'))
    #replace .fcs with .keyword; append .keyword if not ending in .fcs
    if (filepieces[length(filepieces)] == 'fcs') {
      filepieces[length(filepieces)] = 'keyword'
    } else {
      filepieces[length(filepieces)+1] = 'keyword'
    }
    keyword_file <- paste(filepieces, collapse = '.')
  }

  if (debug) {
    print (paste("Converting file: ", input_file))
    print (paste("Original file name: ", filename))
    print (paste("Output file name: ", output_file))
    print (paste("Keyword file name: ", keyword_file))
  }
  fcs <- read.FCS(input_file, transformation=F)
  keywords <- keyword(fcs)
  write.table(as.matrix(keywords), file=keyword_file, quote=F,
              row.names=T, col.names=F, sep='=', append=F)
}

# Extract Keywords
# @param input_file     FCS file to be transformed
# @param keyword_file   FCS file keywords ".keywords" extension"
transformFCS <- function(input_file, keyword_file="", debug=FALSE) {
  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid <- isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })

  if (isValid) {
    extractKeywords(input_file, keyword_file, debug)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}

args <- commandArgs(TRUE)
transformFCS(args[1], args[2])
