#!/usr/bin/env R
VERSION = "0.1"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressPackageStartupMessages(require(RaceID))
suppressPackageStartupMessages(require(FateID))
source(args[1])

do.trajectoryinspection.stemID <- function(ltr){
    makeBranchLink <- function(i,j,k){
        ingoing <- paste(sort(c(i,j)), collapse=".")
        outgoing <- paste(sort(c(j,k)), collapse=".")
        messed <- sort(c(ingoing,outgoing))
        return(list(messed[[1]], messed[[2]]))
    }

    zzz <- do.call(getproj, c(ltr, trjsid.getproj))
    bra <- branchcells(
        ltr,
        do.call("makeBranchLink", as.list(trjsid.branchcells.ijk))
    )
    write.table(
        head(bra$diffgenes$z, trjsid.numdiffgenes),
        file=out.diffgenes)

    print(do.call(plotmap, c(bra$scl, final=TRUE, fr=FALSE)))
    print(do.call(plotmap, c(bra$scl, final=TRUE, fr=TRUE)))
    print(do.call(plotmap, c(bra$scl, final=FALSE, fr=FALSE)))
    print(do.call(plotmap, c(bra$scl, final=FALSE, fr=TRUE)))
}

do.trajectoryinspection.fateID <- function(ltr){
    n <- do.call(cellsfromtree, c(ltr, trjfid.cellsfrom))
    x <- getfdata(ltr@sc)

    trjfid.filterset$x = x
    trjfid.filterset$n = n$f
    fs <- do.call(filterset, c(trjfid.filterset))
    trjfid.getsom$x = fs
    s1d <- do.call(getsom, c(trjfid.getsom))
    trjfid.procsom$s1d = s1d
    ps <- do.call(procsom, c(trjfid.procsom))

    y    <- ltr@sc@cpart[n$f]
    fcol <- ltr@sc@fcol

    trjfid.plotheat$xpart = y
    trjfid.plotheat$xcol = fcol

    ##Plot average z-score for all modules derived from the SOM:
    trjfid.plotheat$x = ps$nodes.z
    trjfid.plotheat$ypart = unique(ps$nodes)
    print(do.call(plotheatmap, c(trjfid.plotheat)))
    ##Plot z-score profile of each gene ordered by SOM modules:
    trjfid.plotheat$x = ps$all.z
    trjfid.plotheat$ypart = ps$nodes
    print(do.call(plotheatmap, c(trjfid.plotheat)))
    ##Plot normalized expression profile of each gene ordered by SOM modules:
    trjfid.plotheat$x = ps$all.e
    trjfid.plotheat$ypart = ps$nodes
    print(do.call(plotheatmap, c(trjfid.plotheat)))
    ##Plot binarized expression profile of each gene (z-score < -1, -1 < z-score < 1, z-score > 1):
    trjfid.plotheat$x = ps$all.b
    trjfid.plotheat$ypart = ps$nodes
    print(do.call(plotheatmap, c(trjfid.plotheat)))

    ## This should be written out, and passed back into the tool
    ## to perform sominspect
    return(list(fs=fs,ps=ps,y=y,fcol=fcol,nf=n$f))
}

do.trajectoryinspection.fateID.sominspect <- function(domo){
    g <- trjfidsomi.use.genes
    if (class(g) == "numeric"){
        g <- names(ps$nodes)[ps$nodes %in% g]
    }

    typ = NULL
    if (!is.null(trjfidsomi.use.types)){
        typ = sub(trjfidsomi.use.types,"", domo$nf)
    }

    trjfidsomi$x = domo$fs
    trjfidsomi$y = domo$y
    trjfidsomi$g = g
    trjfidsomi$n = domo$nf
    trjfidsomi$col = domo$fcol
    trjfidsomi$types = typ

    ## The average pseudo-temporal expression profile of this group
    ## can be plotted by the function plotexpression:
    print(do.call(plotexpression, c(trjfidsomi)))
}

ltr <- in.rdat

pdf(out.pdf)
if (perform.stemID) do.trajectoryinspection.stemID(ltr)
if (perform.fateID) {
    domo <- do.trajectoryinspection.fateID(ltr)
    if (perform.fateID.sominspect) {
        do.trajectoryinspection.fateID.sominspect(domo)
    }
}
dev.off()
