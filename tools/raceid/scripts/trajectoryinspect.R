#!/usr/bin/env R
VERSION <- "0.2" # nolint

args <- commandArgs(trailingOnly = T)

if (length(args) != 1) {
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
suppressWarnings(suppressPackageStartupMessages(require(FateID)))
source(args[1])

test <- list()
test$side <- 3
test$line <- 2.5
second <- test
second$cex <- 0.5
second$line <- 2.5

                                        # nolint start
do.trajectoryinspection.stemID <- function(ltr) {
    makeBranchLink <- function(i, j, k) {
        ingoing <- paste(sort(c(i, j)), collapse=".")
        outgoing <- paste(sort(c(j, k)), collapse=".")
                                        # nolint end
        messed <- sort(c(ingoing, outgoing))
        return(list(messed[[1]], messed[[2]]))
    }

    zzz <- do.call(getproj, c(ltr, trjsid.getproj)) # nolint
    bra <- branchcells(
        ltr,
        do.call("makeBranchLink", as.list(trjsid.branchcells.ijk))
    )
    write.table(
        head(bra$diffgenes$z, trjsid.numdiffgenes),
        file = out.diffgenes)

    par(mfrow = c(2, 2), cex = 0.5)
    print(do.call(plotmap, c(bra$scl, final = FALSE, fr = FALSE))) # nolint
    print(do.call(mtext, c("Initial Clusters (tSNE)", test)))
    print(do.call(plotmap, c(bra$scl, final = TRUE, fr = FALSE))) # nolint
    print(do.call(mtext, c("Final Clusters (tSNE)", test)))
    print(do.call(plotmap, c(bra$scl, final = FALSE, fr = TRUE))) # nolint
    print(do.call(mtext, c("Initial Clusters (F-R)", test)))
    print(do.call(plotmap, c(bra$scl, final = TRUE, fr = TRUE))) # nolint
    print(do.call(mtext, c("Final Clusters (F-R)", test)))
}

do.trajectoryinspection.fateID <- function(ltr) { # nolint
    n <- do.call(cellsfromtree, c(ltr, trjfid.cellsfrom)) # nolint
    x <- getfdata(ltr@sc) # nolint

    trjfid.filterset$x <- x
    trjfid.filterset$n <- n$f # nolint
    fs <- do.call(filterset, c(trjfid.filterset)) # nolint
    trjfid.getsom$x <- fs
    s1d <- do.call(getsom, c(trjfid.getsom)) # nolint
    trjfid.procsom$s1d <- s1d
    ps <- do.call(procsom, c(trjfid.procsom)) # nolint

    y    <- ltr@sc@cpart[n$f] # nolint
    fcol <- ltr@sc@fcol

    trjfid.plotheat$xpart <- y
    trjfid.plotheat$xcol <- fcol

    ##Plot average z-score for all modules derived from the SOM:
    trjfid.plotheat$x <- ps$nodes.z # nolint
    trjfid.plotheat$ypart <- unique(ps$nodes) # nolint
    print(do.call(plotheatmap, c(trjfid.plotheat))) # nolint
    print(do.call(mtext, c(paste(c("Average z-score for all modules ",
                                   "derived from SOM")), test)))
    ##Plot z-score profile of each gene ordered by SOM modules:
    trjfid.plotheat$x <- ps$all.z # nolint
    trjfid.plotheat$ypart <- ps$nodes # nolint
    print(do.call(plotheatmap, c(trjfid.plotheat))) # nolint
    print(do.call(mtext, c(paste(c("z-score profile of each gene ",
                                   "ordered by SOM modules")),
                           test)))
    ## Plot normalized expression profile of each gene
    ## ordered by SOM modules:
    trjfid.plotheat$x <- ps$all.e # nolint
    trjfid.plotheat$ypart <- ps$nodes # nolint
    print(do.call(plotheatmap, c(trjfid.plotheat))) # nolint
    print(do.call(mtext, c(paste(c("Normalized expression profile of each ",
                                   "gene ordered by SOM modules")), test)))
    ## Plot binarized expression profile of each gene
    ## (z-score < -1, -1 < z-score < 1, z-score > 1):
    trjfid.plotheat$x <- ps$all.b # nolint
    trjfid.plotheat$ypart <- ps$nodes # nolint
    print(do.call(plotheatmap, c(trjfid.plotheat))) # nolint
    print(do.call(mtext, c("Binarized expression profile of each gene", test)))
    ## This should be written out, and passed back into the tool
    ## to perform sominspect
    return(list(fs = fs, ps = ps, y = y, fcol = fcol, nf = n$f)) # nolint
}

do.trajectoryinspection.fateID.sominspect <- function(domo) { # nolint
    g <- trjfidsomi.use.genes # nolint
    if (class(g) == "numeric") {
        g <- names(ps$nodes)[ps$nodes %in% g] # nolint
    }

    typ <- NULL
    if (!is.null(trjfidsomi.use.types)) {
        typ <- sub(trjfidsomi.use.types, "", domo$nf) # nolint
    }

    trjfidsomi$x <- domo$fs # nolint
    trjfidsomi$y <- domo$y # nolint
    trjfidsomi$g <- g
    trjfidsomi$n <- domo$nf # nolint
    trjfidsomi$col <- domo$fcol # nolint
    trjfidsomi$types <- typ

    ## The average pseudo-temporal expression profile of this group
    ## can be plotted by the function plotexpression:
    par(mfrow = c(1, 1))
    test$cex <- 1
    second$line <- 1.5
                                        # nolint start
    if (trjfidsomi$name == "Title") trjfidsomi$name <- ""
    print(do.call(plotexpression, c(trjfidsomi)))
    mess2 <- paste(c(trjfidsomi.use.genes), collapse = ", ")
                                        # nolint end
    mess1 <- "Average pseudo-temporal expression profile"
    print(do.call(mtext, c(mess1, test)))
    print(do.call(mtext, c(mess2, second)))
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
