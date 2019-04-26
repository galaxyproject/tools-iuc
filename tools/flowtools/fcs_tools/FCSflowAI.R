#!/usr/bin/Rscript --vanilla
#
# Authors: Gianni Monaco
#
# Reference: flowAI: automatic and interactive anomaly discerning
#            tools for flow cytometry data.
#            Gianni Monaco, Hao Chen, Michael Poidinger, Jinmiao Chen,
#            Joao Pedro de Magalhaes and Anis Larbi
#            Bioinformatics (2016)
#            doi: 10.1093/bioinformatics/btw191
#

library(flowAI)
library(methods)

# parse arguments

args <- commandArgs(trailingOnly = TRUE)

remFS <- if(args[4]) c("FSC", "SSC") else NULL

flow_auto_qc(
    fcsfiles = args[1],
    remove_from = args[2],
    alphaFR = as.numeric(args[3]),
    ChRemoveFS = remFS,
    outlierFS = args[5],
    pen_valueFS = as.numeric(args[6]),
    sideFM = args[7],
    fcs_QC = ifelse(args[9] == "None", FALSE, "_QC"),
    fcs_highQ = ifelse(args[10] == "None", FALSE, "_highQ"),
    fcs_lowQ = ifelse(args[11] == "None", FALSE, "_lowQ"),
    folder_results = FALSE)

try(file.rename(dir(".", pattern = ".*_QC.html"), args[8]), silent =TRUE)
try(file.rename(dir(".", pattern = ".*_QC.fcs"), args[9]), silent =TRUE)
try(file.rename(dir(".", pattern = ".*_highQ.fcs"), args[10]), silent =TRUE)
try(file.rename(dir(".", pattern = ".*_lowQ.fcs"), args[11]), silent =TRUE)
