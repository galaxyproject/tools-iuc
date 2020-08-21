#!/usr/bin/Rscript
# Stacked 1D Density Plot Module for Galaxy
# flowviz
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#

library(flowViz)
library(methods)

generateStackedPlots <- function(fs, chans=list(), output="", flag_pdf=FALSE) {
  h <- 800
  w <- 1200
  if (length(fs@colnames)>8){
    h <- 1200
    w <- 1600
  }
  channels_to_plot <- fs@colnames
  if (length(chans) > 0){
    channels_to_plot <- fs@colnames[chans]
  }

  if (flag_pdf) {
    pdf(output, useDingbats=FALSE, onefile=TRUE)
    print({
      densityplot(~., fs, channels = channels_to_plot)
    })
    dev.off()
  } else {
    png(output, type="cairo", height=h, width=w)
    print({
      densityplot(~., fs, channels = channels_to_plot)
    })
    dev.off()
  }
}

checkFlowSet <- function(fcsfiles, newnames, channels=list(), out_file ="",
                         flag_pdf=FALSE) {

  isValid <- F
  markerCheck <- T

  for (i in 1:length(fcsfiles)){
   is_file_valid <- F
   tryCatch({
     fcs <- read.FCS(fcsfiles[i], transformation=FALSE)
     is_file_valid <- T
   }, error = function(ex) {
     print(paste(ex))
   })
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
  }
  if (markerCheck) {
   isValid <- T
  } else {
   quit(save = "no", status = 12, runLast = FALSE)
  }

  if (isValid) {
    fs <- read.flowSet(files=fcsfiles, transformation=FALSE)
    fs@phenoData@data$name <- newnames
    generateStackedPlots(fs, channels, out_file, flag_pdf)
  } else {
    quit(save = "no", status = 12, runLast = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
channels <- list()
flag_pdf <- FALSE

if (args[1]=="None") {
  flag_default <- TRUE
} else {
  if (args[1] == "i.e.:1,3,4"){
    flag_default <- TRUE
  } else {
    channels <- as.numeric(strsplit(args[1], ",")[[1]])
    for (channel in channels){
      if (is.na(channel)){
        quit(save = "no", status = 11, runLast = FALSE)
      }
    }
  }
}

if (args[3] == "PDF"){
  flag_pdf <- TRUE
}

nb_files <- (length(args) - 3) / 2
fcsfiles <- character(nb_files)
newnames <- character(nb_files)
j <- 1
## get files and file names
for (i in 4:length(args)) {
  if (!i%%2){
    fcsfiles[[j]] <- args[i]
    newnames[[j]] <- args[i+1]
    j <- j + 1
  }
}

checkFlowSet(fcsfiles, newnames, channels, args[2], flag_pdf)
