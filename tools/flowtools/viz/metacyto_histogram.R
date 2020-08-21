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
    quit(save = "no", status = 12, runLast = FALSE)
  } else {
    tmp <- gsub(" ", "", cl_def, fixed = TRUE)
    clean_def <- gsub(",", "|", tmp, fixed = TRUE)
    return (clean_def)
  }
}

generatePlots <- function(fpath="", fname="", gates=vector(), outdir="", uc="",
                          flag_pdf=F){
  dir.create(outdir)
  ff <- read.FCS(fpath, truncate_max_range=F)
  markers <- markerFinder(ff)
  colnames(ff@exprs) <- markers

  sc <- searchCluster(fcsFrame=ff, clusterLabel=gates)

  if (length(gates) == length(sc$clusterList)){
    sink(uc)
    cat("All provided cluster definition were used.")
    sink()
  } else {
    unused_cluster <- setdiff(gates, names(sc$clusterList))
    write.table(unused_cluster, uc, quote=F, row.names=F, col.names = F)
  }

  groupname <- unlist(strsplit(fname, ".fcs"))[[1]]
  extension <- if (flag_pdf) "plot.pdf" else "plot.png"
  for (i in 1:length(sc$clusterList)){
    gate <- gsub("|","", names(sc$clusterList[i]), fixed=T)
    plotname <- paste(c(groupname, gate, extension), collapse="_")
    outplot <- file.path(outdir, plotname)
    if (flag_pdf){
      pdf(outplot, useDingbats=F, onefile=T)
      par(mfrow=c(2,2))
      for (j in 1:length(markers)) {
        if (markers[[j]]!="SAMPLE_ID" && markers[[j]]!="TIME") {
          plot_title <- paste0(markers[[j]],", cluster definition:\n", gate)
          x_all <- ff@exprs[,markers[[j]]]
          b <- seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
          subset <- ff@exprs[sc$clusterList[[i]], markers[[j]]]
          hist(x_all, col=rgb(0,0,0), xlab=markers[[j]], breaks=b, freq=T,
               border=F, main=plot_title)
          hist(subset, add=T, breaks=b, col=rgb(1,0,0), freq=T, border=F)
          if (markers[[j]] %in% names(sc$cutoff)){
            abline(v=sc$cutoff[markers[[j]]])
          }
        }
      }
      dev.off()
    } else {
      markers_ct <- length(markers) - length(grep(x=markers, pattern="SAMPLE_ID|TIME"))
      nb_rows <- ceiling(markers_ct / 2)
      h <- nb_rows*400
      png(outplot, type="cairo", height=h, width=800)
      par(mfrow=c(nb_rows,2))
      for (j in 1:length(markers)) {
        if (markers[[j]]!="SAMPLE_ID" && markers[[j]]!="TIME") {
          plot_title <- paste0(markers[[j]],", cluster definition:\n", gate)
          x_all <- ff@exprs[,markers[[j]]]
          b <- seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
          subset <- ff@exprs[sc$clusterList[[i]],markers[[j]]]
          hist(x_all, col=rgb(0,0,0), xlab=markers[[j]], breaks=b, freq=T,
               border=F, main=plot_title)
          hist(subset, add=T, breaks=b, col=rgb(1,0,0), freq=T, border=F)
          if (markers[[j]] %in% names(sc$cutoff)){
            abline(v=sc$cutoff[markers[[j]]])
          }
        }
      }
      dev.off()
    }
  }
}

checkFCSfile <- function(inputf="", inputn="", clusters=vector(),
                         output_dir="", unused="", flag=F){
  isValid <- FALSE
  tryCatch({
    isValid <- isFCSfile(inputf)
  }, error = function(ex) {
    print(paste("Input file is not a valid FCS file.", ex))
  })
  if (isValid) {
    generatePlots(inputf, inputn, clusters, output_dir, unused, flag)
  } else {
    quit(save = "no", status = 12, runLast = FALSE)
  }
}

################################################################################
################################################################################
args <- commandArgs(trailingOnly = TRUE)

gates <- vector()
if (args[6] == "F"){
  ## obvs deal with it if file
  cluster_file <- read.table(args[7], header=F,colClasses="character")
  gates <- unlist(cluster_file)
} else {
  cl_df <- args[7:length(args)]
  gates <- sapply(cl_df, checkClusterDef)
}

flag_pdf <- if (args[5]=="PDF") TRUE else FALSE
gate_list <- toupper(gates)
checkFCSfile(args[1], args[2], gate_list, args[3], args[4], flag_pdf)
