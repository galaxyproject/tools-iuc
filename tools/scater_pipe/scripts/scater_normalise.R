#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

## Load libs, common functions, source Galaxy config
source(paste(script_dir, "common.R", sep="/"))

## Already have sce here
## message("method"); message(nmethod)

sce <- normaliseExprs(sce, method = nmethod, return_log = ret_log, return_norm_as_exprs = T)
saveRDS(sce, "sce_out.rds")
