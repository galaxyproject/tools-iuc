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

library(MetaCyto)

checkPanel <- function(df, outfile="", pdf_out="") {
  report <- panelSummary(df,".",cluster=F,width=30,height=20)
  if (pdf_out!=""){
    file.rename("./panel_summary.pdf", pdf_out)
  }
  markers <- data.frame('Markers'=row.names(report))
  s <- cbind(markers, report)
  write.table(s, file=outfile, quote=F, row.names=F, col.names=T, sep='\t')
}

checkInputFormat <- function(infile="", outfile="", pdf_file="") {
  df <- read.table(infile, sep="\t", header=T, colClasses="character")
  nm <- colnames(df)
  check_ab <- if ("antibodies" %in% nm) TRUE else FALSE
  check_sdy <- if ("study_id" %in% nm) TRUE else FALSE

  if (check_sdy && check_ab){
    checkPanel(df, outfile, pdf_file)
  } else {
    quit(save = "no", status = 10, runLast = FALSE)
  }
}

################################################################################
args <- commandArgs(trailingOnly = TRUE)

pdf_out <- ""
if (length(args)==3){
  pdf_out <- args[3]
}

checkInputFormat(args[1], args[2], pdf_out)
