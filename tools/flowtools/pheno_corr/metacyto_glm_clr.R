#!/usr/bin/env Rscript
######################################################################
#                  Copyright (c) 2018 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1 - April 2018
# Author: Cristel Thomas
#
#

library(MetaCyto)
library(data.table)
library(dplyr)
library(tidyr)

maxNb <- function(cluster_defs) {
  max <- 0
  for (l in cluster_defs) {
    m <- length(strsplit(l, "\\|")[[1]])
    if (m>max) {max <- m}
  }
  return(max)
}

collectDataTSV <- function(files=vector()) {
  all_stats <- NULL
  for (fn in files){
    cs <- read.table(fn, sep="\t", header=T, stringsAsFactors=FALSE)
    nms <- colnames(cs)
    ## this is going to break at some point:
    w <- which(nms=="fcs_names") + 1
    t1 <- w:ncol(cs)
    cluster_stat <- gather(cs, parameter_name, value, !!t1)
    all_stats <- rbind(all_stats, cluster_stat)
  }
  return(all_stats)
}

checkClusterStat <- function(cs_file){
  cs <- read.table(cs_file, sep="\t", header=T, colClasses="character")
  nm <- toupper(colnames(cs))
  to_check <- c("FRACTION", "GROUP_NAME", "LABEL")
  fp <- c("FCS_FILES", "FCS_NAMES")
  t1 <- sum(to_check %in% nm)
  t2 <- sum(fp %in% nm)
  if (t1 == 3 && t2 > 1){
    colnames(cs) <- nm
    filenames <- if ("FCS_NAMES" %in% nm) unique(data.frame(cs$GROUP_NAME, cs$FCS_NAMES)) else unique(data.frame(cs$GROUP_NAME, cs$FCS_FILES))
    return (filenames)
  } else {
    quit(save = "no", status = 11, runLast = FALSE)
  }
}

findCluster <- function(cs_file, gate){
  cs <- read.table(cs_file, sep="\t", header=T, colClasses="character")
  if (gate %in% cs$label){
    return(1)
  } else {
    return(0)
  }
}

runMA <- function(cs_files=vector(), md=list(), params=vector(), output_ma="",
                   outplot="", ci=0.95, re="", loo="", clr=""){

  fcs_stats <- collectDataTSV(cs_files)
  all_data <- inner_join(fcs_stats, md, by=c("fcs_names", "group_name"))
  other_vars <- if (length(params)==1) NULL else params[2:length(params)]

  dat=subset(all_data,all_data$parameter_name=="fraction"&
             all_data$label==clr)

  if (outplot != "") {
    pdf(outplot, useDingbats=F, onefile=T, height=7, width=7)
  }
  MA <- metaAnalysis(value="value",
                  variableOfInterst=params[[1]],
                  main=clr,
                  otherVariables=other_vars,
                  studyID="group_name",
                  data=dat,
                  CILevel=ci,
                  ifScale=c(T,F))

  if (outplot != "") {
    dev.off()
  }

  colnames(MA$Summary)[[1]] <- "group_name"
  write.table(MA$Summary, file=output_ma, quote=F, row.names=F, col.names=T, sep="\t")
  write.table(MA$leave1out, file=loo, quote=F, row.names=F, col.names=T, sep="\t")
  sink(re)
  print(MA$REmodel)
  sink()
}


checkInput <- function(stat_files=vector(), md_file="", params=NULL,
                       output_ma="", outplot="", ci=0.95, re="", loo="",
                       clr_def="", flag_default=F){

  allfiles <- rbindlist(lapply(stat_files, checkClusterStat))
  colnames(allfiles) <- c("group", "filename")

  md <- read.table(md_file, sep="\t", header=T, colClasses="character")
  md_names <- toupper(colnames(md))
  gp_check <- c("GROUP_NAME", "STUDY_ID")
  t3 <- sum(gp_check %in% md_names)
  if (t3 != 1) {
    quit(save = "no", status = 13, runLast = FALSE)
  }
  colnames(md) <- md_names

  idx_files <- grep("FCS_FILES", md_names)
  if (length(idx_files) == 0) {
    idx_files <- grep("FCS_NAMES", md_names)
    if (length(idx_files)==0) {
      idx_files <- grep("FILENAMES", md_names)
      if (length(idx_files)==0) {
        quit(save = "no", status = 12, runLast = FALSE)
      }
    }
  }
  idx_gp <- if ("GROUP_NAME" %in% md_names) grep("GROUP_NAME", md_names) else grep("STUDY_ID", md_names)
  md_fn <- unique(data.frame(md[,idx_gp], md[,idx_files]))

  # one line per file:group
  if (dim(md_fn)[[1]] != dim(md)[[1]]) {
    quit(save = "no", status = 14, runLast = FALSE)
  }
  colnames(md_fn) <- c("group", "filename")
  colnames(md)[idx_files] <- "fcs_names"
  colnames(md)[idx_gp] <- "group_name"
  allfiles_test <- paste(allfiles$group, allfiles$filename, sep="##")
  md_fn_test <- paste(md_fn$group, md_fn$filename, sep="##")
  if (length(intersect(allfiles_test, md_fn_test)) != length(allfiles_test)) {
    quit(save = "no", status = 15, runLast = FALSE)
  }

  main_param <- NULL
  if (flag_default){
    main_param <- colnames(md)[[3]]
  } else {
    if (!is.null(params)){
      t4 <- sum(params %in% md_names)
      if (t4 != length(params)){
        quit(save = "no", status = 16, runLast = FALSE)
      }
    }
  }

  all_params <- c(main_param, params)

  # check that cluster definition exists in the list of clusters
  # if not exit 17
  # if in only one group exit 18
  cluster_check <- sapply(stat_files, findCluster, gate=clr_def)
  if (sum(cluster_check) == 0) {
    quit(save = "no", status = 17, runLast = FALSE)
  } else if (sum(cluster_check) == 1){
    quit(save = "no", status = 18, runLast = FALSE)
  }

  runMA(stat_files, md, all_params, output_ma, outplot, ci, re, loo, clr_def)
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)

flag_default <- if (args[2]==""||args[2]=="i.e.: Gender"||args[2]=="None") T else F
ci <- as.numeric(args[4])

cluster <- ""
if (args[5]=="" || args[5]=="None" || args[5]=="i.e.:CD3+,CD4-,CD8+,CCR7+"){
  quit(save = "no", status = 19, runLast = FALSE)
} else {
  tmp <- gsub(" ", "", args[5], fixed = TRUE)
  clean_def <- gsub(",", "|", tmp, fixed = TRUE)
  cluster <- toupper(clean_def)
}

i <- grep(args, pattern="STAT_FILES")

outplot <- ""
extra_params <- NULL

if (i==10){
  outplot <- args[9]
} else if (i>8){
  if (args[8]=="" || args[8]=="i.e.: Age, Ethnicity" || args[8]=="None") {
    quit(save = "no", status = 10, runLast = FALSE)
  }
  extras <- gsub(", ", ",", args[8], fixed = TRUE)
  extra_params <- strsplit(toupper(extras), ",")[[1]]
  if (i==11) {
    outplot <- args[10]
  }
}

params <- if (flag_default) extra_params else c(toupper(args[2]), extra_params)

ii <- i + 1
stat_files <- args[ii:length(args)]

checkInput(stat_files, args[1], params, args[3], outplot, ci, args[6], args[7],
           cluster, flag_default)
