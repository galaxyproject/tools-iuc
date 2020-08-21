#!/usr/bin/env Rscript
######################################################################
#                  Copyright (c) 2018 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1 - January 2018
# Author: Cristel Thomas
#
#

library(flowCore)
library(MetaCyto)

compareLists <- function(m1, m2){
  list_check <- T
  if (is.na(all(m1==m2))) {
    mm1 <- is.na(m1)
    mm2 <- is.na(m2)
    if (all(mm1==mm2)){
      if (!all(m1==m2, na.rm=TRUE)){
        listCheck <- F
      }
    } else {
      list_check <- F
    }
  } else if (!all(m1==m2)) {
    list_check <- F
  }
  return(list_check)
}


runBatchProcessing <- function(sampling_size=5000, flag_default=T,
                               to_exclude, outdir="", outfile="",
                               labels, assays, factors, fcspaths, fcsnames) {
  # Create meta_data object
  fp <- unlist(fcspaths)
  file_counts <- lengths(fcspaths)
  group_names <- rep(labels,times=file_counts)
  group_bs <- rep(factors,times=file_counts)
  group_types <- rep(assays,times=file_counts)

  meta_data <- data.frame("fcs_files"=fp, "study_id"=group_names)

  # excluded_parameters
  default_param <- c("FSC-A", "FSC-H", "FSC-W", "FSC", "SSC-A", "SSC-H",
                     "SSC-W", "SSC", "Time", "Cell_length", "cell_length",
                     "CELL_LENGTH")
  excluded_parameters <- if (flag_default) default_param else to_exclude
  # Run preprocessing.batch
  preprocessing.batch(inputMeta=meta_data,
                      assay=group_types,
                      b=group_bs,
                      fileSampleSize=sampling_size,
                      outpath=outdir,
                      excludeTransformParameters = excluded_parameters)

  # deal with outputs
  # output[2]: a csv file summarizing the pre-processing result.
  ## -> open file to pull info and print out filenames rather than path.
  tmp_csv <- file.path(outdir, "processed_sample_summary.csv")
  tmp <- read.csv(tmp_csv)
  tmp$old_index <- seq(1, length(tmp$fcs_names))

  fn <- unlist(fcsnames)
  df <- data.frame("fcs_files"=fp, "filenames"=fn)

  # merge two data frames by ID
  total <- merge(tmp,df,by="fcs_files")
  total2 <- total[order(total$old_index),]
  to_drop <- c("fcs_names", "fcs_files", "old_index")
  newdf <- total2[ , !(names(total2) %in% to_drop)]
  write.table(newdf, file=outfile, quote=F, row.names=F, col.names=T, sep='\t')

  file.remove(tmp_csv)
}

checkFCS <- function(sampling=5000, flag_default=TRUE, to_exclude,
                     outdir="", outfile="", labels, assays, factors,
                     fcspaths, fcsnames) {

  marker_pb <- FALSE
  for (i in 1:length(fcspaths)){
    for (n in 1:length(fcspaths[[i]])) {
      marker_check <- FALSE
      marker_channel <- FALSE
      is_file_valid <- FALSE
      tryCatch({
        fcs <- read.FCS(fcspaths[[i]][[n]], transformation=FALSE)
        is_file_valid <- TRUE
      }, error = function(ex) {
        print(paste("File is not a valid FCS file:", fnames[[i]][[n]] , ex))
      })

      if (is_file_valid){
        if (n == 1) {
          m1 <- as.vector(pData(parameters(fcs))$desc)
          c1 <- colnames(fcs)
        } else {
          m2 <- as.vector(pData(parameters(fcs))$desc)
          c2 <- colnames(fcs)
          marker_check <- compareLists(m1,m2)
          marker_channel <- compareLists(c1,c2)
        }
        if (n>1 && marker_check==F) {
          marker_pb <- TRUE
          print(paste("Marker discrepancy detected in markers -- group", labels[[i]]))
        } else if (n>1 && marker_channel==F) {
          marker_pb <- TRUE
          print(paste("Marker discrepancy detected in channels -- group", labels[[i]]))
        }
      } else {
        quit(save = "no", status = 10, runLast = FALSE)
      }
    }
  }

  if (marker_pb){
    quit(save = "no", status = 12, runLast = FALSE)
  } else {
    runBatchProcessing(sampling, flag_default, to_exclude, outdir, outfile,
                       labels, assays, factors, fcspaths, fcsnames)
  }
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)

sub_sampling <- NULL
if (as.numeric(args[1]) > 0) {
  sub_sampling <- as.numeric(args[1])
}

# parameters to exclude => args[4]
to_exclude <- vector()
flag_default <- FALSE
i <- 1
if (args[4]=="None" || args[4]== "" || args[4] == "i.e.:FSC, SSC, CD88") {
  flag_default <- TRUE
} else {
  excluded <- strsplit(args[4], ",")[[1]]
  for (channel in excluded){
    stripped_chan <- gsub(" ", "", channel, fixed = TRUE)
    if (is.na(stripped_chan)==FALSE){
      to_exclude[[i]] <- stripped_chan
    }
    i <- i + 1
  }
}

tot_args <- length(args)
tmpargs <- paste(args[5:tot_args], collapse="=%=")
tmppanels <- strsplit(tmpargs, "=%=new_panel=%=")
nb_panel <- length(tmppanels[[1]])

labels <- vector(mode="character", length=nb_panel)
assay_types <- vector(mode="character", length=nb_panel)
scaling_factors <- vector(mode="numeric", length=nb_panel)
filepaths <- list()
filenames <- list()

j <- 1
for (pnl in tmppanels[[1]]) {
  tmppanel <- strsplit(pnl, "=%=")
  nb_files <- (length(tmppanel[[1]]) - 3 ) /2
  tmplist <- character(nb_files)
  tmpnames <- character(nb_files)
  if (tmppanel[[1]][[1]] == "None" || tmppanel[[1]][[1]] == ""){
    quit(save = "no", status = 11, runLast = FALSE)
  } else {
    labels[[j]] <- tmppanel[[1]][[1]]
  }
  # assay type
  assay_types[[j]] <- tmppanel[[1]][[2]]

  scaling_factors[[j]] <- 0
  if (as.numeric(tmppanel[[1]][[3]]) > 0) {
    scaling_factors[[j]] <- 1 / as.numeric(tmppanel[[1]][[3]])
  }

  k <- 1
  for (m in 4:length(tmppanel[[1]])){
    if (!m%%2){
      tmplist[[k]] <- tmppanel[[1]][[m]]
      tmpnames[[k]]<- tmppanel[[1]][[m+1]]
      k <- k + 1
    }
  }
  filepaths[[tmppanel[[1]][1]]] <- tmplist
  filenames[[tmppanel[[1]][1]]] <- tmpnames
  j <- j + 1
}

checkFCS(sub_sampling, flag_default, to_exclude, args[2], args[3], labels,
         assay_types, scaling_factors, filepaths, filenames)
