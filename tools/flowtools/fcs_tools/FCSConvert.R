#!/usr/bin/Rscript --vanilla
# ImmPort FCSConvert
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Converts the FCS file to text without transformaton
# To run in R
# 1) library(flowCore)
# 2) source("FCSConvert.R")
# 3) transformFCS("filename")
#
# Version 1.4.1
# March 2016 -- added lines to run directly from command line
#

library(flowCore)

convertFCS <- function(fcs,compensate=FALSE,debug=FALSE) {
  # Check file type and FCS version
  if (class(fcs)[1] != "flowFrame") {
    print("convertFCS requires flowFrame object as input")
    return(FALSE)
  }

  keywords <- keyword(fcs)
  markers <- colnames(fcs)
  print_markers <- as.vector(pData(parameters(fcs))$desc)
  # Update print_markers if the $P?S not in the FCS file
  for (i in 1:length(print_markers)) {
    if (is.na(print_markers[i])) {
      print_markers[i] <- markers[i]
    }
  }

  if (debug) {
    print("****Inside convertFCS")
    print(paste("    FCS version:", keywords$FCSversion))
    print(paste("    DATATYPE:", keywords['$DATATYPE']))
  }

  if (keywords$FCSversion == "2" ||
      keywords$FCSversion == "3" ||
      keywords$FCSversion == "3.1" ) {
      datatype = unlist(keywords['$DATATYPE'])
    if (datatype == 'F') {
        # Apply compensation if available and requested
      spill <- keyword(fcs)$SPILL
      if (is.null(spill) == FALSE && compensate == TRUE) {
        if (debug) {
          print("Attempting compensation")
        }
        tryCatch({ fcs = compensate(fcs, spill) },
                   error = function(ex) { str(ex); })
      }
      # Process fcs expression data, using transformation
      # based on category of the marker.
      fcs_exprs <- exprs(fcs)
      colnames(fcs_exprs) <- print_markers
        } else if (datatype == 'I') {
          fcs_exprs <- exprs(fcs)
          colnames(fcs_exprs) = print_markers
        } else {
          print(paste("Data type", datatype, "in FCS 3 is not supported"))
          fcs_exprs <- FALSE
        }
  } else {
    print(paste("FCS version", keyword(fcs)$FCSversion, "is not supported"))
    fcs_exprs <- FALSE
  }
  fcs_exprs
}


#
# Starting function for processing a FCS file
#
processFCSFile <- function(input_file, output_file="",
                           keyword_file="",compensate=FALSE, debug=FALSE) {

  #
  # Generate the file names for the output_file and keyword_file
  #
  pieces <- unlist(strsplit(input_file, .Platform$file.sep))
  filename <- pieces[length(pieces)]

  if (output_file == "") {
    filepieces = unlist(strsplit(filename, '\\.'))
    #replace .fcs with .txt; append .txt if not ending in .fcs
    if (filepieces[length(filepieces)] == 'fcs') {
      filepieces[length(filepieces)] <- 'txt'
    } else {
      filepieces[length(filepieces)+1] <- 'txt'
    }
    output_file <- paste(filepieces, collapse = '.')
  }

  if (keyword_file == "") {
    filepieces <- unlist(strsplit(filename, '\\.'))
    #replace .fcs with .keyword; append .keyword if not ending in .fcs
    if (filepieces[length(filepieces)] == 'fcs') {
      filepieces[length(filepieces)] <- 'keyword'
    } else {
      filepieces[length(filepieces)+1] <- 'keyword'
    }
    keyword_file <- paste(filepieces, collapse = '.')
  }

  if (debug) {
    print (paste("Converting file: ",input_file))
    print (paste("Original file name: ",filename))
    print (paste("Output file name: ",output_file))
    print (paste("Keyword file name: ",keyword_file))
  }
  fcs <- read.FCS(input_file, transformation=F)
  keywords <- keyword(fcs)
  write.table(as.matrix(keywords),file=keyword_file, quote=F,
              row.names=T, col.names=F, sep='=', append=F)

  transformed_data <- convertFCS(fcs,compensate,debug)
  write.table(transformed_data, file=output_file, quote=F,
              row.names=F,col.names=T, sep='\t', append=F)
}

# Convert FCS file without transformation
# @param input_file     FCS file to be transformed
# @param output_file    FCS file transformed ".txt" extension
# @param keyword_file   FCS file keywords ".keywords" extension"
# @param compensate     Flag indicating whether to apply compensation
#                       matrix if it exists.
transformFCS <- function(input_file, output_file="",
                         compensate=FALSE, keyword_file="", debug=FALSE) {

  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid <- isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })

  if (isValid) {
    processFCSFile(input_file, output_file, keyword_file, compensate, debug)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
transformFCS(args[1], args[2], args[3])
