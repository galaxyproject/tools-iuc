#!/usr/bin/Rscript --vanilla
######################################
#  Support for FCS sniffer function  #
######################################

library(flowCore)

checkFCSfile <- function(inputf) {
  isValid <- FALSE
  tryCatch({
    isValid <- isFCSfile(inputf)
  }, error = function(ex) {
    print("error reading file")
  })
  if (isValid) {
    print(TRUE)
  } else {
    print("not a valid FCS file")
  }
}
args <- commandArgs(trailingOnly = TRUE)

checkFCSfile(args[1])
