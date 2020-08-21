#!/usr/bin/Rscript
# Aggregate FCS files Module for Galaxy
# with FlowSOM AggregateFlowFrames
######################################################################
#                  Copyright (c) 2017 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#
library(FlowSOM)
library(flowCore)

downsampleMergeFCS <- function(fcs_files, nb_cells, output="", flag_fcs=FALSE) {
  ff <- AggregateFlowFrames(fcs_files, nb_cells, writeOutput = FALSE)
  n <- length(colnames(ff)) - 2
  exprs(ff) <- exprs(ff)[,1:n]
  if (flag_fcs) {
    write.FCS(ff, output)
  } else {
    saveRDS(ff, file=output)
  }
}

checkFCSfiles <- function(fcsfiles, ds_factor=0.1, out_file ="",
                         flag_fcs=FALSE) {
  isValid <- F
  nb_events <- 0
  markerCheck <- T

  for (i in 1:length(fcsfiles)){
    is_file_valid <- F
    tryCatch({
      fcs <- read.FCS(fcsfiles[i], transformation=FALSE)
      is_file_valid <- T
      nb_events <- nb_events + as.numeric(fcs@description$`$TOT`)
    }, error = function(ex) {
      print(paste(ex))
    })
    if (is_file_valid){
      if (i == 1) {
        m1 <- as.vector(pData(parameters(fcs))$desc)
      } else {
        m2 <- as.vector(pData(parameters(fcs))$desc)
        if (is.na(all(m1==m2))) {
          mm1 <- is.na(m1)
          mm2 <- is.na(m2)
          if (all(mm1==mm2)){
            if (!all(m1==m2, na.rm=TRUE)){
              markerCheck <- F
            }
          } else {
            markerCheck <- F
          }
        } else if (!all(m1==m2)) {
          markerCheck <- F
        }
      }
    } else {
      quit(save = "no", status = 10, runLast = FALSE)
    }
  }

  if (markerCheck) {
    isValid <- T
  } else {
    quit(save = "no", status = 12, runLast = FALSE)
  }

  ## translate ds_factor to nb of events
  nb_cell <- floor(ds_factor*nb_events)

  if (isValid) {
    downsampleMergeFCS(fcsfiles, nb_cell, out_file, flag_fcs)
  } else {
    quit(save = "no", status = 10, runLast = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
flag_fcs <- FALSE

if (args[2] == "FCS"){
  flag_fcs <- TRUE
}

if (args[3] == "" || args[3] == "i.e.:0.1 or 10X") {
  factor <- 0.1
} else {
  #rm last X if it's there
  ds <- gsub("X", "", args[3])
  if (!is.na(as.numeric(ds))) {
    factor <- as.numeric(ds)
    if (factor > 1 && factor <= 100) {
      factor <- as.numeric(ds) / 100
    } else if (factor > 100){
      quit(save = "no", status = 11, runLast = FALSE)
    }
  } else {
    quit(save = "no", status = 11, runLast = FALSE)
  }
}

nb_files <- (length(args) - 3)
fcsfiles <- character(nb_files)
j <- 1
## get files and file names
for (i in 4:length(args)) {
    fcsfiles[[j]] <- args[i]
    j <- j + 1
}

checkFCSfiles(fcsfiles, factor, args[1], flag_fcs)
