#!/usr/bin/Rscript --vanilla
# FCS Headers Module for Galaxy
# FlowCore
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 2
# May 2018
# Cristel Thomas
#
#

library(flowCore)

getFCSChannels <- function(input_fcs) {
  fcs <- read.FCS(input_fcs, transformation=F)
  return(colnames(fcs))
}

getFCSMarkers <- function(input_fcs){
  ffcs <- read.FCS(input_fcs, transformation=F)
  fmarkers <- as.vector(pData(parameters(ffcs))$desc)
  return(fmarkers)
}

getFCSMarkerNames <- function(output_file="", file_paths=vector(),
                              fcs_names=vector()) {
  check_files <- sapply(file_paths, isFCSfile)
  channels <- lapply(file_paths[check_files], getFCSChannels)
  markers <- lapply(file_paths[check_files], getFCSMarkers)

  nb_col <- max(lengths(channels))
  nc <- lapply(channels, `length<-`, nb_col)
  ct <- t(as.data.frame(nc))

  nm <- lapply(markers, `length<-`, nb_col)
  mt <- t(as.data.frame(nm))

  nb_files <- sum(check_files)
  Index <- rep(c("channels", "markers"), each=nb_files)
  Filename <- rep(fcs_names[check_files], 2)

  idx_nb <- seq(nb_col)
  ttt <- rbind(ct,mt)
  finalt <- cbind(Filename, Index, ttt)
  colnames(finalt)[3:length(colnames(finalt))] <- idx_nb


  if (nb_files != length(file_paths)){
    not_fcs <- fcs_names[!check_files]
    new_df <- cbind(not_fcs, "Not a valid FCS file")
    empty_frame <- data.frame(matrix("", nrow=length(not_fcs), ncol=nb_col),
                              stringsAsFactors = F)
    not_fcs_files <- cbind(new_df, empty_frame)
    colnames(not_fcs_files) <- colnames(finalt)
    new_final <- rbind(finalt, not_fcs_files)
    write.table(new_final, file=output_file, quote=F, row.names=F, col.names=T,
                sep='\t', append=F)
    quit(save = "no", status = 10, runLast = FALSE)
  } else {
    # output file
    write.table(finalt, file=output_file, quote=F, row.names=F, col.names=T,
                sep='\t',
    append=F)
  }
}

################################################################################
################################################################################

args <- commandArgs(trailingOnly = TRUE)

i <- 1
nb_files <- (length(args)-1) / 2
fcs_files <- character(nb_files)
fcs_names <- character(nb_files)
for (j in 1:length(args)){
  if (!j%%2){
    fcs_files[[i]] <- args[[j]]
    fcs_names[[i]]<- args[[j+1]]
    i <- i + 1
  }
}

getFCSMarkerNames(args[1], fcs_files, fcs_names)
