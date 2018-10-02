#!/usr/bin/env R
VERSION = "0.1"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
source(args[1])

do.pseudotemp <- function(sc){
    ltr <- Ltree(sc)
    ltr <- compentropy(ltr)
    ltr <- do.call(projcells, c(ltr, pstc.projc))
    ltr <- do.call(projback, c(ltr, pstc.projb))
    ltr <- lineagegraph(ltr)
    ltr <- do.call(comppvalue, c(ltr, pstc.comppval))
    print(do.call(plotgraph, c(ltr, pstc.plotgraph)))
    x <- do.call(compscore, c(ltr, pstc.compscore))
    ## make sure this gets plotted
    ##Plots
    pdf(out.pdf)
    print(plotdistanceratio(ltr))
    print(plotspantree(ltr))
    print(plotprojections(ltr))
    print(plotlinkscore(ltr))
    print(projenrichment(ltr))
    dev.off()
    return(ltr)
}

ltr <- do.pseudotemp(in.rdat)
saveRDS(ltr, out.rdat)
