#!/usr/bin/Rscript
# Module for Galaxy
# Compares groups of FCS to FlowSOM reference tree
# with FlowSOM
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

checkFiles <- function(groups){
  all_files <- unlist(groups)
  all_unique <- unique(all_files)
  if (length(all_unique) != length(all_files)) {
      quit(save = "no", status = 14, runLast = FALSE)
  }
}

compareLists <- function(m1, m2){
  listCheck <- T
  if (is.na(all(m1==m2))) {
    mm1 <- is.na(m1)
    mm2 <- is.na(m2)
    if (all(mm1==mm2)){
      if (!all(m1==m2, na.rm=TRUE)){
        listCheck <- F
      }
    } else {
      listCheck <- F
    }
  } else if (!all(m1==m2)) {
    listCheck <- F
  }
  return(listCheck)
}

prettyMarkerNames <- function(flowFrame){
  n <- flowFrame@parameters@data[, "name"]
  d <- flowFrame@parameters@data[, "desc"]
  d[is.na(d)] <- n[is.na(d)]
  prettyNames <- list()
  if(any(grepl("#",d))){
      # Support for hashtag notation:
      # antibody#fluorochrome -> antibody (fluorochrome)
      prettyNames <- gsub("#(.*)$"," (\\1)",d)
  } else {
      prettyNames <- paste(d, " <", n, ">", sep="")
  }
  return(prettyNames)
}

compareToTree <- function(fst, wilc_thresh=0.05, output="", plot="", stats,
                          comp_groups, filenames) {
  groupRes <- CountGroups(fst, groups=comp_groups, plot=FALSE)
  pdf(plot, useDingbats=FALSE, onefile=TRUE)
  tresh <- wilc_thresh
  pg <- PlotGroups(fst, groupRes, p_tresh=tresh)
  dev.off()

  nb_nodes <- length(pg[[1]])
  nb_comp <- length(pg)
  m <- matrix(0, nrow=nb_nodes, ncol=nb_comp+1)
  s <- seq_len(nb_nodes)
  m[,1]<- as.character(s)
  for (i in 1:nb_comp){
    m[s,i+1]<- as.character(pg[[i]])
  }
  groupnames <- attr(comp_groups,"names")
  out_colnames <- paste(groupnames, collapse="-")
  colnames(m) <- c("Node",out_colnames)
  write.table(m, file=output, quote=F, row.names=F, col.names=T, sep='\t',
              append=F)

  ## get filenames
  filepaths <- unlist(comp_groups)
  fnames <- unlist(filenames)
  nb_files <- length(filepaths)
  comp_files <- list()
  for (i in 1:length(filepaths)){
    comp_files[[filepaths[[i]]]] <- fnames[[i]]
  }

  group_list <- list()
  for (grp in attr(comp_groups, "names")) {
    for (f in comp_groups[[grp]]){
      group_list[[f]] <- grp
    }
  }
  out_stats <- attr(stats, "names")
  if ("counts" %in% out_stats){
    gp_counts <- as.matrix(groupRes$counts)
    tpc <- matrix("", nrow=nb_files, ncol=2)
    tpc[,1] <- as.character(lapply(rownames(gp_counts), function(x) comp_files[[x]]))
    tpc[,2] <- as.character(lapply(rownames(gp_counts), function(x) group_list[[x]]))
    gp_counts <- cbind(tpc, gp_counts)
    colnames(gp_counts)[[1]] <- "Filename"
    colnames(gp_counts)[[2]] <- "Group"
    t_gp_counts <- t(gp_counts)
    write.table(t_gp_counts, file=stats[["counts"]], quote=F, row.names=T, col.names=F, sep='\t',
                append=F)
  }
  if ("pctgs" %in% out_stats){
    gp_prop <- as.matrix(groupRes$pctgs)
    tpp <- matrix("", nrow=nb_files, ncol=2)
    tpp[,1] <- as.character(lapply(rownames(gp_prop), function(x) comp_files[[x]]))
    tpp[,2] <- as.character(lapply(rownames(gp_prop), function(x) group_list[[x]]))
    gp_prop <- cbind(tpp, gp_prop)
    colnames(gp_prop)[[1]] <- "Filename"
    colnames(gp_prop)[[2]] <- "Group"
    t_gp_prop <- t(gp_prop)
    write.table(t_gp_prop, file=stats[["pctgs"]], quote=F, row.names=T, col.names=F, sep='\t',
                append=F)
  }
  if ("means" %in% out_stats){
    gp_mean <- as.matrix(groupRes$means)
    t_gp_mean <- t(gp_mean)
    tpm <- matrix(0, nrow=nb_nodes, ncol=1)
    tpm[,1] <- seq_len(nb_nodes)
    t_gp_mean <- cbind(tpm, t_gp_mean)
    colnames(t_gp_mean)[[1]] <- "Nodes"
    write.table(t_gp_mean, file=stats[["means"]], quote=F, row.names=F, col.names=T, sep='\t',
                append=F)
  }
  if ("medians" %in% out_stats){
    gp_med <- as.matrix(groupRes$medians)
    t_gp_med <- t(gp_med)
    tpd <- matrix(0, nrow=nb_nodes, ncol=1)
    tpd[,1] <- seq_len(nb_nodes)
    t_gp_med <- cbind(tpd, t_gp_med)
    colnames(t_gp_med)[[1]] <- "Nodes"
    write.table(t_gp_med, file=stats[["medians"]], quote=F, row.names=F, col.names=T, sep='\t',
                append=F)
  }
}

checkFCS <- function(tree, output="", plot="", thresh = 0.05, stats, groups,
                     filenames) {

  fcsfiles <- unlist(groups)
  tree_valid <- F
  markerCheck <- T
  tryCatch({
    fsomtree <- readRDS(tree)
    tree_valid <- T
  }, error = function(ex) {
    print(paste(ex))
  })

  fst <- if (length(fsomtree)==2) fsomtree[[1]] else fsomtree

  if (tree_valid){
    tree_markers <- as.vector(fst$prettyColnames)
    tree_channels <- as.vector(colnames(fst$data))
    if (length(tree_markers) < 1){
      quit(save = "no", status = 11, runLast = FALSE)
    }
  } else {
    quit(save = "no", status = 11, runLast = FALSE)
  }

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
      c1 <- colnames(fcs)
      # compare to tree markers
      pm <- prettyMarkerNames(fcs)
      if (!all(tree_markers %in% pm)){
        quit(save = "no", status = 13, runLast = FALSE)
      }
    } else {
      m2 <- as.vector(pData(parameters(fcs))$desc)
      c2 <- colnames(fcs)
      markerCheck <- compareLists(m1,m2)
      markerChannel <- compareLists(c1,c2)
    }
  }
  if (markerCheck && markerChannel) {
    compareToTree(fst, thresh, output, plot, stats, groups, filenames)
  } else {
    quit(save = "no", status = 12, runLast = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)

first_g1 <- 5
tot_args <- length(args)
g <- list()
tmplist <- c("counts", "means", "medians", "pctgs")

for (i in 5:13){
  if (args[i] %in% tmplist){
    first_g1 <- first_g1 + 2
    g[[args[i]]] <- args[i+1]
  }
}

tmpargs <- paste(args[first_g1:tot_args], collapse="=%=")
tmpgroups <- strsplit(tmpargs, "=%=DONE=%=")

groups <- list()
filenames <- list()
for (gps in tmpgroups[[1]]) {
  tmpgroup <- strsplit(gps, "=%=")
  nb_files <- (length(tmpgroup[[1]]) - 1 ) /2
  tmplist <- character(nb_files)
  tmpnames <- character(nb_files)
  j <- 1
  for (i in 2:length(tmpgroup[[1]])){
    if (!i%%2){
      tmplist[[j]] <- tmpgroup[[1]][i]
      tmpnames[[j]]<- tmpgroup[[1]][i+1]
      j <- j + 1
    }
  }
  groups[[tmpgroup[[1]][1]]] <- tmplist
  filenames[[tmpgroup[[1]][1]]] <- tmpnames
}

checkFiles(groups)
checkFCS(args[1], args[2], args[3], args[4], g, groups, filenames)
