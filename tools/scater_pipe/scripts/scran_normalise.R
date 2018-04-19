#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

## Load libs, common functions, source Galaxy config
source(paste(script_dir, "common.R", sep="/"))

## Already have sce here
qclust <- NULL

if (!is.null(qclust_minsize)){
    qclust <- quickCluster(sce, qclust_minsize)
}

sce_new <- computeSumFactors(sce, clusters = qclust, positive = c_positive, min.mean = c_min_mean)
sce_scran <- normalise(sce)

saveRDS(sce_scran, "sce_out.rds")
