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

generateTree <- function(ff, output="", columns=list(), cluster=10, xgrid=10,
                         ygrid=10,plot="", plot_pdf=FALSE, mplot="", flag_def=T,
                         table="", mtable="", flag_meta=FALSE, user_seed=42,
                         flag_nodesize=F) {

  # check default -- if def get all except FSC/SSC
  # also check nb of markers/channels and indices
  markers <- colnames(ff)
  print_markers <- as.vector(pData(parameters(ff))$desc)
  # Update print_markers if the $P?S not in the FCS file
  for (i in 1:length(print_markers)) {
    if (is.na(print_markers[i])) {
      print_markers[i] <- markers[i]
    }
  }

  if (flag_def){
    channels_to_exclude <- c(grep(markers, pattern="FSC"),
                             grep(markers, pattern="SSC"),
                             grep(markers, pattern="Time"))
    columns <- markers[-channels_to_exclude]
  }

  set.seed(user_seed)
  fs <- ReadInput(ff, compensate=F, transform=F, scale=T)
  fs <- BuildSOM(fs, colsToUse = columns, xdim=xgrid, ydim=ygrid)
  fst <- BuildMST(fs, tSNE=T)

  if (!mplot==""){
    pdf(mplot, useDingbats=FALSE, onefile=TRUE)
    for (marker in markers){
      PlotMarker(fst, marker)
    }
    dev.off()
  }
  metaC <- metaClustering_consensus(fst$map$codes, k=cluster, seed=user_seed)

  if (!plot==""){
    if (flag_nodesize){
      fst <- UpdateNodeSize(fst, reset=TRUE)
      fst$MST$size <- fst$MST$size/2
    }
    if (plot_pdf) {
      pdf(plot, useDingbats=FALSE)
      PlotStars(fst, backgroundValues = as.factor(metaC))
      dev.off()
    } else {
      png(plot, type="cairo", height=800, width=800)
      PlotStars(fst, backgroundValues = as.factor(metaC))
      dev.off()
    }
  }
  if (!table==""){
    m <- matrix(0,nrow=nrow(ff),ncol=1)
    s <- seq_len(nrow(ff))
    if (flag_meta){
      m[s,] <- metaC[fst$map$mapping[,1]]
    } else {
      m[s,] <- fst$map$mapping[,1]
    }
    colnames(m) <- "FlowSOM"
    ff <- cbind2(ff,m)
    out <- exprs(ff)
    print_markers <- append(print_markers, "Population")
    colnames(out) <- print_markers
    write.table(out, file=table, quote=F, row.names=F, col.names=T, sep='\t',
                append=F)

    nb_nodes <- max(fst$map$mapping[,1])
    mm <- matrix(0, nrow=nb_nodes, ncol=2)
    ss <- seq_len(nb_nodes)
    mm[,1]<- as.character(ss)
    mm[ss,2]<- as.character(metaC)
    colnames(mm) <- c("Node", "Meta-Cluster")
    write.table(mm, file=mtable, quote=F, row.names=F, col.names=T, sep='\t',
                append=F)

  }
  saveRDS(fst, file = output)
}

flowFrameOrFCS <- function(input, output="", columns=list(),cluster=10,xgrid=10,
                           ygrid=10,plot="",plot_pdf=FALSE, mplot="", default=T,
                           table="", mtable="", flag_meta=FALSE, user_seed=42,
                           nodesize=FALSE) {
  isValid <- F
  is_fcs <- F
  is_ff <- F
  ff <- ""
  tryCatch({
    is_fcs <- isFCSfile(input)
  }, error = function(ex) {
    print(paste(ex))
  })

  if (!is_fcs){
    tryCatch({
      ff <- readRDS(input)
      is_ff <- T
    }, error = function(ex) {
      print(paste(ex))
    })
  } else {
    ff <- read.FCS(input, transformation=FALSE)
  }

  if (!is_ff && !is_fcs) {
    quit(save = "no", status = 10, runLast = FALSE)
  } else {
    for (cols in columns){
      if (cols > length(colnames(ff))){
        quit(save = "no", status = 12, runLast = FALSE)
      }
    }
    generateTree(ff, output, columns, cluster, xgrid, ygrid, plot, plot_pdf,
                 mplot, default, table, mtable, flag_meta, user_seed, nodesize)
  }
}

args <- commandArgs(trailingOnly = TRUE)
flag_default <- FALSE
columns <- list()

if (args[3] == "" || args[3] == "i.e.:1,2,5") {
  flag_default <- TRUE
} else {
  #rm last X if it's there
  columns <- as.numeric(strsplit(args[3], ",")[[1]])
  for (col in columns){
    if (is.na(col)){
      quit(save = "no", status = 11, runLast = FALSE)
    }
  }
}

cluster <- 10
if (!args[4] == ""){
  if (!is.na(as.integer(args[4]))){
    cluster <- as.integer(args[4])
  } else {
    quit(save = "no", status = 13, runLast = FALSE)
  }
}

xgrid <- 10
if (!args[5] == ""){
  if (!is.na(as.integer(args[5]))){
    cluster <- as.integer(args[5])
  } else {
    quit(save = "no", status = 14, runLast = FALSE)
  }
}

ygrid <- 10
if (!args[6] == ""){
  if (!is.na(as.integer(args[6]))){
    cluster <- as.integer(args[6])
  } else {
    quit(save = "no", status = 14, runLast = FALSE)
  }
}
seed <- 42
if (!args[7]==""){
  if (!is.na(as.integer(args[7]))){
    seed <- as.integer(args[7])
  } else {
    quit(save = "no", status = 15, runLast = FALSE)
  }
}

plot <- ""
mplot <- ""
plot_pdf <- FALSE
table <- ""
mtable <- ""
flag_meta <- FALSE
nodesize <- FALSE
nb_args <- length(args)

if (nb_args==16) {
  plot <- args[8]
  if (args[9]=='PDF') {
    plot_pdf <- TRUE
  }
  nodesize <- args[10]
  mplot <- args[11]
  table <- args[13]
  mtable <- args[14]
  if (args[12]=='meta'){
    flag_meta<-TRUE
  }
} else if (nb_args==15){
  plot <- args[8]
  if (args[9]=='PDF') {
    plot_pdf <- TRUE
  }
  nodesize <- args[10]
  table <- args[12]
  mtable <- args[13]
  if (args[11]=='meta'){
    flag_meta<-TRUE
  }
} else if (nb_args==13) {
  mplot <- args[8]
  table <- args[10]
  mtable <- args[11]
  if (args[9]=='meta'){
    flag_meta<-TRUE
  }
} else if (nb_args==12) {
  table <- args[9]
  mtable <- args[10]
  if (args[8]=='meta'){
    flag_meta<-TRUE
  }
} else if (nb_args==11) {
  plot <- args[8]
  if (args[9]=='PDF') {
    plot_pdf <- TRUE
  }
  nodesize <- args[10]
  mplot <- args[11]
} else if (nb_args==10) {
  plot <- args[8]
  if (args[9]=='PDF') {
    plot_pdf <- TRUE
  }
  nodesize <- args[10]
} else if (nb_args==8){
  mplot <- args[8]
}

flowFrameOrFCS(args[1], args[2], columns, cluster, xgrid, ygrid, plot, plot_pdf,
               mplot, flag_default, table, mtable, flag_meta, seed, nodesize)
