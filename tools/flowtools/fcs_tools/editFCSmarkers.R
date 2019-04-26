#!/usr/bin/Rscript
# modify channels and marker names in FCS
#
######################################################################
#                  Copyright (c) 2017 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Cristel Thomas
# Version 2 - May 2018
# Modified to take in marker/channel names by name rather than index
#
library(flowCore)

checkCandM <- function(m_set, channels=vector(), markers=vector()){
  if (m_set[[3]]){
    in_file <- m_set[[1]] %in% channels
  } else {
    in_file <- m_set[[1]] %in% markers
  }
  yay <- if (sum(in_file)==0) F else T
  return(yay)
}

modifyMarkersFCS <- function(input, output="", report="", flag_fcs=F,
                             marker_sets=list()) {

  fcs <- read.FCS(input, transformation=F)
  original_channels <- colnames(fcs)
  original_markers <- as.vector(pData(parameters(fcs))$desc)
  nb <- length(original_channels)
  ## check if markers are in FCS files
  check_markers <- sapply(marker_sets, checkCandM, channels=original_channels,
                          markers=original_markers)
  if (sum(check_markers)==0) {
    quit(save = "no", status = 13, runLast = FALSE)
  }

  post_channels <- colnames(fcs)
  post_markers <- as.vector(pData(parameters(fcs))$desc)

  for (m_set in marker_sets) {
    if (m_set[[3]]){
      chan_to_replace <- post_channels %in% m_set[[1]]
      for (i in 1:nb){
        if (chan_to_replace[[i]]){
          post_channels[[i]] <- m_set[[2]]
        }
      }
    } else {
      marker_to_replace <- post_markers %in% m_set[[1]]
      for (i in 1:nb){
        if (marker_to_replace[[i]]){
          post_markers[[i]] <- m_set[[2]]
          pm <- paste("$P", as.character(i), "S", sep="")
          fcs@description[[pm]] <- m_set[[2]]          
        }
      }
    }
  }

  colnames(fcs) <- post_channels
  pData(parameters(fcs))$desc <- post_markers

  # write report
  sink(report)
  cat("###########################\n")
  cat("##    BEFORE RENAMING    ##\n")
  cat("###########################\nFCS Channels\n")
  cat("---------------------------\n")
  cat(original_channels,"---------------------------", "FCS Markers","---------------------------",original_markers, sep="\n")
  cat("\n###########################\n")
  cat("##    AFTER  RENAMING    ##\n")
  cat("###########################\nFCS Channels\n")
  cat("---------------------------\n")
  cat(post_channels,"---------------------------","FCS Markers","---------------------------", post_markers, sep="\n")
  sink()

  # output fcs
  if (flag_fcs) {
    write.FCS(fcs, output)
  } else {
    saveRDS(fcs, file = output)
  }
}

checkFCS <- function(fcsfile, out_file ="", report="", flag_fcs=FALSE,
                     marker_sets=list()) {
  isValid <- F
  tryCatch({
    isValid <- isFCSfile(fcsfile)
  }, error = function(ex) {
    print(paste(ex))
  })

  if (isValid) {
    modifyMarkersFCS(fcsfile, out_file, report, flag_fcs, marker_sets)
  } else {
    quit(save = "no", status = 10, runLast = FALSE)
  }
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)
flag_fcs <- if (args[3]=="FCS") TRUE else FALSE

items <- args[5:length(args)]
marker_sets <- list()
j <- 1
for (i in seq(1, length(items), 3)) {
  if (items[i]=="None" || items[i]== "" || items[i]== "i.e.:TLR 6, TLR6PE") {
    quit(save = "no", status = 11, runLast = FALSE)
  }
  if (items[i+1]=="None" || items[i+1]=="" || items[i+1]=="i.e.:TLR6") {
    quit(save = "no", status = 12, runLast = FALSE)
  }

  old_names <- strsplit(items[i], ",")[[1]]
  to_replace <- sapply(old_names, function(x) trimws(x))
  replacement <- trimws(items[i+1])
  flag_channel <- if (items[i+2]=="C") TRUE else FALSE
  m_set <- list(to_replace, replacement, flag_channel)
  marker_sets[[j]] <- m_set
  j <- j + 1
}

checkFCS(args[1], args[2], args[4], flag_fcs, marker_sets)
