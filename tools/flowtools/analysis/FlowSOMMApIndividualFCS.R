#!/usr/bin/Rscript
# Module for Galaxy
# Generates FlowSOM reference tree
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

## geometric mean from
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

prettyMarkerNames <- function(flowFrame){
  n <- flowFrame@parameters@data[, "name"]
  d <- flowFrame@parameters@data[, "desc"]
  d[is.na(d)] <- n[is.na(d)]
  prettyNames <-list()
  if(any(grepl("#",d))){
      # Support for hashtag notation:
      # antibody#fluorochrome -> antibody (fluorochrome)
      prettyNames <- gsub("#(.*)$"," (\\1)",d)
  } else {
      prettyNames <- paste(d, " <", n, ">", sep="")
  }
  return(prettyNames)
}

checkMarkers <- function(fcsfiles, flag_ff=FALSE){
  markerCheck <- T
  for (i in 1:length(fcsfiles)){
    if (flag_ff){
      fcs <- readRDS(fcsfiles[i])
    } else {
      fcs <- read.FCS(fcsfiles[i], transformation=FALSE)
    }
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
  if (!markerCheck) {
    quit(save = "no", status = 13, runLast = FALSE)
  }
}

mapToTree <- function(filenames, filepaths, flag_ff=FALSE, reftree, cluster=10,
                      outdir="",flag_meta=FALSE, mfi='mfi', stat1="", stat2="",
                      stat3="", plot="", mplot="") {

  checkMarkers(filepaths, flag_ff)

  # get tree
  fst <- readRDS(reftree)
  plots <- FALSE
  mplots <- FALSE
  dir.create(outdir)
  if (!plot==""){
    dir.create(plot)
    plots <- TRUE
  }
  if (!mplot==""){
    dir.create(mplot)
    mplots <- TRUE
  }
  metaC <- metaClustering_consensus(fst$map$codes, k=cluster, seed=33)
  nb_pop <- if (flag_meta) cluster else max(fst$map$mapping[,1])
  nb_samples <- length(filepaths)
  nb_marker <- length(fst$prettyColnames)
  print_markers <- gsub(' <.*>','',fst$prettyColnames)
  print_markers_ff <- append(print_markers, "Population")

  m_stat1 <- matrix(0, nrow=nb_samples, ncol=(nb_pop+2))
  colnames(m_stat1) <- c("FileID", "SampleName", seq_len(nb_pop))

  sink(stat2)
  cat(print_markers, sep="\t")
  cat("\tPercentage\tPopulation\tSampleName\n")
  sink()

  col_m3 <- c("Population")
  for (m in print_markers){
    m1 <- paste(m, "mean", sep="_")
    m2 <- paste(m, "median", sep="_")
    m3 <- paste(m, "stdev", sep="_")
    col_m3 <- append(col_m3, c(m1,m2,m3))
  }
  col_stat3 <- c(col_m3, "Percentage_mean","Percentage_median","Percentage_stdev")

  for (i in 1:length(filepaths)){
    if (flag_ff){
      ff <- readRDS(filepaths[i])
    } else {
      ff <- read.FCS(filepaths[i], transformation=FALSE)
    }
    if (i==1) {
      # compare to tree markers
      pm <- prettyMarkerNames(ff)
      if (!all(fst$prettyColnames %in% pm)){
        quit(save = "no", status = 14, runLast = FALSE)
      }
    }

    fsom <- NewData(fst, ff)

    if (mplots){
      markers <- colnames(ff)
      tmpmplot <- paste(filenames[i], "marker_plots.pdf", sep="_")
      pdf(file.path(mplot,tmpmplot), useDingbats=FALSE, onefile=TRUE)
      for (marker in markers){
        PlotMarker(fsom, marker)
      }
      dev.off()
    }

    if (!plot==""){
      plotpath <- paste(filenames[i], "tree.png", sep="_")
      png(file.path(plot, plotpath), type="cairo", height=800, width=800)
      PlotStars(fsom, backgroundValues = as.factor(metaC))
      dev.off()
    }

    m <- matrix(0,nrow=nrow(ff),ncol=1)
    s <- seq_len(nrow(ff))
    if (flag_meta){
      m[s,] <- metaC[fsom$map$mapping[,1]]
    } else {
      m[s,] <- fsom$map$mapping[,1]
    }
    colnames(m) <- "FlowSOM"
    ff <- cbind2(ff,m)
    out <- exprs(ff)
    colnames(out) <- print_markers_ff

    clr_table <- paste(filenames[i], "clustered.flowclr", sep="_")
    write.table(out, file=file.path(outdir, clr_table), quote=F, row.names=F, col.names=T, sep='\t',
                append=F)

    cluster_count <- table(out[,"Population"])
    cluster_prop <- prop.table(cluster_count)*100
    m1_tmp <- numeric(nb_pop)
    for (j in 1:nb_pop){
      if (as.character(j) %in% names(cluster_count)){
        m1_tmp[[j]] <- format(round(cluster_prop[[as.character(j)]], 2), nsmall=2)
      }
    }
    samplename <- paste("Sample", i, sep="")
    m_stat1[i,] <- c(filenames[i], samplename, m1_tmp)
    # flowstat2
    # Marker1 Marker2 Marker3 ... Population Percentage SampleName
    # MFIs for each marker
    # dimension ==> col = nb of markers + 3; row = nb of files * nb of clusters
    if (mfi=="mfi"){
      m2 <- aggregate(out[,1:nb_marker], list(out[,nb_marker+1]), mean)
    } else if (mfi=="mdfi") {
      m2 <- aggregate(out[,1:nb_marker], list(out[,nb_marker+1]), median)
    } else {
      m2 <- aggregate(out[,1:nb_marker], list(out[,nb_marker+1]), gm_mean)
    }

    m2["Percentage"] <- as.character(cluster_prop)
    m2["Population"] <- as.character(m2$Group.1)
    m2["SampleName"] <- samplename
    m2t <- as.matrix(m2[2:length(m2)])
    write.table(m2t, file=stat2, quote=F, row.names=F, col.names=F, sep='\t',
                append=T)

  }
  write.table(m_stat1, file=stat1, quote=F, row.names=F, col.names=T, sep='\t',
              append=F)

  m2df <- read.table(stat2, sep="\t", header=T)
  ag <- aggregate(m2df[,0:nb_marker+1], list(m2df[,nb_marker+2]), function(x) c(mn=mean(x), med=median(x), stdv=sd(x)))
  m3t <- as.matrix(ag)
  colnames(m3t) <- col_stat3
  write.table(m3t, file=stat3, quote=F, row.names=F, col.names=T, sep='\t',
              append=F)
}

flowFrameOrFCS <- function(filenames, filepaths, reftree, cluster=10, outdir="",
                           flag_meta=FALSE, mfi='mfi', stat1="", stat2="",
                           stat3="", plot="", mplot=""){

  isValid <- F
  is_fcs <- F
  is_ff <- F
  flag_ff <- FALSE
  i <- 0
  for (f in filepaths){
    tryCatch({
      is_fcs <- isFCSfile(f)
    }, error = function(ex) {
      print(paste(ex))
    })
    if (!is_fcs){
      tryCatch({
        ff <- readRDS(f)
        is_ff <- T
      }, error = function(ex) {
        print(paste(ex))
      })
    } else {
      i <- i + 1
    }
    if (!is_ff && !is_fcs) {
      quit(save = "no", status = 10, runLast = FALSE)
    }
  }
  if (i==0) {
    flag_ff <- TRUE
  } else if (!i==length(filenames)){
    quit(save = "no", status = 12, runLast = FALSE)
  }
  mapToTree(filenames, filepaths, flag_ff, reftree, cluster, outdir, flag_meta,
            mfi, stat1, stat2, stat3, plot, mplot)
}

args <- commandArgs(trailingOnly = TRUE)
plot <- ""
mplot <- ""
m <- 8
flag_meta <- FALSE
if (args[4]=='meta') {
  flag_meta <- TRUE
}
if (args[9] == 'newDataTrees') {
  plot <- 'newDataTrees'
  m <- m + 1
  if (args[10] == 'newDataMarkers') {
    mplot <- "newDataMarkers"
    m <- m + 1
  }
} else if (args[9] == 'newDataMarkers') {
  mplot <- "newDataMarkers"
  m <- m + 1
}

n <- m+1
# files + filenames
nb_files <- (length(args) - m) / 2
files1 <- character(nb_files)
files2 <- character(nb_files)
j <- 1
file_list <- args[n:length(args)]
for (i in 1:length(file_list)) {
  if (i%%2){
    files1[[j]] <- file_list[i]
    files2[[j]] <- file_list[i+1]
    j <- j + 1
  }
}

flowFrameOrFCS(files2, files1, args[1], as.numeric(args[3]), args[2], flag_meta,
               args[5], args[6], args[7], args[8], plot, mplot)
