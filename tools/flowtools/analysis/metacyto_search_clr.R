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

checkClusterDef <- function(cl_def){
  if (cl_def=="" || cl_def=="None" || cl_def=="i.e.:CD3+,CD4-,CD8+,CCR7+"){
    quit(save = "no", status = 14, runLast = FALSE)
  } else {
    tmp <- gsub(" ", "", cl_def, fixed = TRUE)
    clean_def <- gsub(",", "|", tmp, fixed = TRUE)
    return (toupper(clean_def))
  }
}

pathToGroupFile <- function(path_to_result) {
  grp <- basename(dirname(path_to_result))
  return(paste(grp, "fcs", sep=".", collapse=NULL))
}

groupFileToGroupName <- function(result_file){
  return(strsplit(result_file, ".", fixed=TRUE)[[1]][1])
}


searchClusterPanels <- function(df, fcspaths, fcsnames, outdir="", uc="",
                              clusters=vector()) {

  working_dir <- "tmp_metacyto"
  working_out <- "tmp_metacyto_out"
  dir.create(working_dir)
  dir.create(outdir)

  # reformat summary -- expects csv + 'fcs_names' && 'fcs_files'
  new_df <- file.path(working_dir, "processed_sample_summary.csv")
  df$fcs_names <- df$filenames
  df$fcs_files <- df$filenames
  write.csv(df, file=new_df, row.names=F)

  # move && rename FCS files to same directory
  for (i in 1:length(fcspaths)) {
    new_file <- file.path(working_dir, fcsnames[[i]])
    file.copy(fcspaths[[i]], new_file)
  }

  searchCluster.batch(preprocessOutputFolder=working_dir,
                      outpath=working_out,
                      clusterLabel=clusters)

  result_files <- list.files(working_out,
                             pattern="cluster_stats_in_each_sample",
                             recursive=T,
                             full.names=T)

  nb_groups <- length(fcsnames)
  no_results <- vector()
  if (length(result_files) != nb_groups) {
   groups_with_results <- sapply(result_files, pathToGroupFile)
   ## one or more groups with no results, figure out which
   no_results <- setdiff(fcsnames, groups_with_results)
  }

  if (length(no_results)==nb_groups){
    sink(uc)
    cat("No clusters were found in none of the groups.")
    sink()
  } else {
    unused_clrs <- list()

    if (length(no_results>0)) {
      grp_no_results <- sapply(no_results, groupFileToGroupName)
      unused_clrs <- data.frame("cluster_label"="any", "not_found_in"=grp_no_results)
    }

    for (result in result_files) {
      group_name <- strsplit(result, .Platform$file.sep)[[1]][2]
      new_filename <- paste(c(group_name, "cluster_stats.txt"), collapse="_")
      new_path <- file.path(outdir, new_filename)
      tmp_df <- read.csv(result)

      used_clr <- as.character(unique(tmp_df$label))
      if (length(used_clr) != length(clusters)) {
        unused <- setdiff(clusters, used_clr)
        tmp_udf <- data.frame("cluster_label"=unused, "not_found_in"=group_name)
        unused_clrs <- rbind(unused_clrs, tmp_udf)
      }
      colnames(tmp_df)[[1]] <- "group_name"
      write.table(tmp_df, new_path, quote=F, row.names=F, col.names=T, sep="\t")
    }

    if (is.null(dim(unused_clrs))){
      sink(uc)
      cat("All provided cluster definition were found in all provided FCS files.")
      sink()
    } else {
      write.table(unused_clrs, uc, quote=F, row.names=F, col.names = T, sep="\t")
    }
  }
}


checkInput <- function(report="", outdir="", list_unused="", list_clusters="",
                       fcs_files=list(), grp_names=list(), clusters=vector()) {
  # check FCS files
  fcspaths <- unlist(fcs_files)
  fcsnames <- unlist(grp_names)
  ct_files <- 0
  some_pb <- FALSE
  for (i in 1:length(fcspaths)){
    is_file_valid <- FALSE
    tryCatch({
      fcs <- read.FCS(fcspaths[[i]], transformation=FALSE)
      is_file_valid <- TRUE
    }, error = function(ex) {
      print(paste("File is not a valid FCS file:", fcsnames[[i]] , ex))
    })
    if (is_file_valid){
      metacyto_pp_check <- if ("sample_id" %in% colnames(fcs)) TRUE else FALSE
      if (metacyto_pp_check) {
        idx <- length(colnames(fcs))
        ct_files <- ct_files + max(fcs@exprs[,idx])
      } else {
        quit(save = "no", status = 11, runLast = FALSE)
      }
    } else {
      some_pb <- TRUE
    }
  }
  # check summary file format
  df <- read.table(report, sep="\t", header=T, colClasses="character")
  nm <- colnames(df)
  check_ab <- if ("antibodies" %in% nm) TRUE else FALSE
  check_sdy <- if ("study_id" %in% nm) TRUE else FALSE

  if (check_sdy && check_ab){
    # check that summary index compatible with FCSs in collection - by number of files == index nb
    if (ct_files != length(df$antibodies)){
      quit(save = "no", status = 12, runLast = FALSE)
    }
  } else {
    quit(save = "no", status = 13, runLast = FALSE)
  }

  if (some_pb){
    quit(save = "no", status = 10, runLast = FALSE)
  } else {
    write.table(clusters, list_clusters, quote=F, row.names=F, col.names = F)
    searchClusterPanels(df, fcspaths, fcsnames, outdir, list_unused, clusters)
  }
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)

i <- grep(args, pattern="FCS_FILES")

cluster_def <- vector()
cl_df <- args[3]
if (i>6){
  ii <- i-1
  more_cl <- args[6:ii]
  cl_df <- c(cl_df, more_cl)
}
cluster_def <- sapply(cl_df, checkClusterDef)

fcs_files <- list()
fcs_names <- list()
j <- 1
m <- i+1
tmp_fcs <- args[m:length(args)]

for (k in 1:length(tmp_fcs)){
  if (k%%2){
    fcs_files[[j]] <- tmp_fcs[[k]]
    fcs_names[[j]]<- tmp_fcs[[k+1]]
    j <- j + 1
  }
}

checkInput(args[1], args[2], args[4], args[5], fcs_files, fcs_names,
           cluster_def)
