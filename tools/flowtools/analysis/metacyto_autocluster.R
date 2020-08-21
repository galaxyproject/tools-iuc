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

autoClusterPanels <- function(params, df, fcspaths, fcsnames, quant=0.95,
                              events=0.05, cluster_algorithm="FlowSOM",
                              clusters=vector(), outdir="", list_clust="",
                              metacluster=40, xdim=10, ydim=10, seed=42, uc=""){

  working_dir <- "tmp_metacyto"
  working_out <- "tmp_metacyto_out"
  dir.create(working_dir)
  dir.create(outdir)

  # get nb of groups
  nb_groups <- length(fcsnames)

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

#### will need to add other parameters when Zicheng has them working.
  if (cluster_algorithm=="FlowSOM") {
    cluster_label <- autoCluster.batch(preprocessOutputFolder=working_dir,
                                        excludeClusterParameters=params,
                                        labelQuantile=quant,
                                        clusterFunction=flowSOM.MC,
                                        minPercent=events,
                                        k=metacluster,
                                        xdim=xdim,
                                        ydim=ydim,
                                        seed=seed)

  } else {
    cluster_label <- autoCluster.batch(preprocessOutputFolder=working_dir,
                                        excludeClusterParameters=params,
                                        labelQuantile=quant,
                                        clusterFunction=flowHC,
                                        minPercent=events)
  }


  # Add potential user-defined label to cluster definitions
  if (length(clusters)>1){
    cluster_label <- c(cluster_label, clusters)
  }
  write.table(cluster_label, list_clust, quote=F, row.names=F, col.names = F)

  # Derive summary statistics for the clusters
  # Result will be written out to the directory speficied by the "outpath" argument
  searchCluster.batch(preprocessOutputFolder=working_dir,
                      outpath=working_out,
                      clusterLabel=cluster_label)

  result_files <- list.files(working_out,
                             pattern="cluster_stats_in_each_sample",
                             recursive=T,
                             full.names=T)
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
      if (length(used_clr) != length(cluster_label)) {
        unused <- setdiff(cluster_label, used_clr)
        tmp_udf <- data.frame("cluster_label"=unused, "not_found_in"=group_name)
        unused_clrs <- rbind(unused_clrs, tmp_udf)
      }
      colnames(tmp_df)[[1]] <- "group_name"
      write.table(tmp_df, new_path, quote=F, row.names=F, col.names=T, sep="\t")
    }

    if (is.null(dim(unused_clrs))){
      sink(uc)
      cat("All provided cluster definition were found in provided FCS files.")
      sink()
    } else {
      write.table(unused_clrs, uc, quote=F, row.names=F, col.names = T, sep="\t")
    }
  }
}

checkInput <- function(params=vector(), report="", fcs_files=list(),
                       grp_names=list(), quant=0.95, events=0.05,
                       cluster_algorithm="FlowSOM",clusters=vector(), outdir="",
                       list_clust="", metacluster=40, xdim=10, ydim=10,
                       seed=42, unused="") {
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
    # check that summary index compatible with FCSs in collection - by number of files == index nb?
    if (ct_files != length(df$antibodies)){
      quit(save = "no", status = 12, runLast = FALSE)
    }
  } else {
    quit(save = "no", status = 13, runLast = FALSE)
  }

  if (some_pb){
    quit(save = "no", status = 10, runLast = FALSE)
  } else {
    autoClusterPanels(params, df, fcspaths, fcsnames, quant, events,
                      cluster_algorithm, clusters, outdir, list_clust,
                      metacluster, xdim, ydim, seed, unused)
  }
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)

ex_param <- c("FSC-A", "FSC-W", "FSC-H", "FSC", "SSC", "SSC-A", "SSC-W",
              "SSC-H", "Time", "Cell_length", "cell_length", "CELL_LENGTH")

if (args[5] != "" && args[5] != "None" && args[5] != "i.e.:FSC,SSC,CD88"){
  tmp <- gsub(" ", "", args[5], fixed = TRUE)
  eparam <- unlist(strsplit(tmp, ","))
  ex_param <- toupper(eparam)
}

i <- grep(args, pattern="PARAM")
ii <- grep(args, pattern="FCS_FILES")

cluster_def <- vector()
if (i>8){
  id <- i-1
  cl_df <- args[8:id]
  cluster_def <- sapply(cl_df, checkClusterDef)
}

metacluster <- 40
xdim <- 10
ydim <- 10
seed <- 42
if (i+1 != ii) {
  metacluster <- as.numeric(args[i+1])
  xdim <- as.numeric(args[i+2])
  ydim <- as.numeric(args[i+3])
  seed <- as.numeric(args[i+4])
}

fcs_files <- list()
fcs_names <- list()
j <- 1
m <- ii+1
n <- length(args) - 1
tmp_fcs <- args[m:n]

for (k in 1:length(tmp_fcs)){
  if (k%%2){
    fcs_files[[j]] <- tmp_fcs[[k]]
    fcs_names[[j]]<- tmp_fcs[[k+1]]
    j <- j + 1
  }
}

checkInput(ex_param, args[1], fcs_files, fcs_names, as.numeric(args[3]),
           as.numeric(args[4]), args[2], cluster_def, args[6], args[7],
           metacluster, xdim, ydim, seed, args[length(args)])
