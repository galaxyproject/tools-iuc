#!/usr/bin/env R
VERSION = "0.2"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
source(args[1])

## layout
test <- list()
test$side = 3
test$line = 3

do.plotting <- function(sc){
    print(plotmap(sc, final = FALSE, fr = FALSE))
    print(do.call(mtext, c("Initial Clustering tSNE", test)))
    print(plotmap(sc, final = TRUE, fr = FALSE))
    print(do.call(mtext, c("Final Clustering tSNE", test)))
    print(plotmap(sc, final = FALSE, fr = TRUE))
    print(do.call(mtext, c("Initial Clustering Fruchterman-Reingold", test)))
    print(plotmap(sc, final = TRUE, fr = TRUE))
    print(do.call(mtext, c("Final Clustering Fruchterman-Reingold", test)))
}


do.inspect.symbolmap <- function(sc){
    if (!is.null(plotsym.use.typeremoveregex)){
        plotsym$types = sub(plotsym.use.typeremoveregex, "", colnames(sc@ndata))

        if (!is.null(plotsym.use.typeremoveregex.subselect)){
            plotsym$subset = plotsym$types[grep(plotsym.use.typeremoveregex.subselect, plotsym$types)]
        }
    }
    plotsym$fr = FALSE
    print(do.call(plotsymbolsmap, c(sc, plotsym)))
    print(do.call(mtext, c("Symbols tSNE", test)))
    plotsym$fr = TRUE
    print(do.call(plotsymbolsmap, c(sc, plotsym)))
    print(do.call(mtext, c("Symbols FR", test)))
}

do.inspect.diffgene <- function(sc){

    getSubNames <- function(lob, sc){
        use.names <- NULL
        if (!is.null(lob$manual)){
            use.names <- lob$manual
        }
        else if (!is.null(lob$regex)){
            nm <- colnames(sc@ndata)
            use.names <- nm[grep(lob$regex, nm)]
        }
        else if (!is.null(lob$cln)){
            use.names <- names(sc@cpart)[sc@cpart %in% lob$cln]
        }
        if (is.null(use.names)){
            stop("A or B names not given!")
        }
        return(use.names)
    }

    A <- getSubNames(gfdat.A.use, sc)
    B <- getSubNames(gfdat.B.use, sc)

    fdat <- getfdata(sc, n=c(A,B))
    dexp <- diffexpnb(fdat, A=A, B=B)
    ## options for diffexpnb are mostly about DESeq, ignore
    plotdiffg$x = dexp
    print(do.call(plotdiffgenesnb, c(plotdiffg)))
    print(do.call(mtext, c("Diff Genes", test)))
}


do.inspect.genesofinterest <- function(sc){
    if (is.null(plotexp$n)){ ## No title, and one gene? Use gene name
        if (length(plotexp$g) == 1){
            plotexp$n <- plotexp$g
        } else {
            plotexp$n <- paste(plotexp$g, collapse=", ")
        }
    }

    title <- paste(":", plotexp$n)
    plotexp$n <- ""

    plotexp$logsc=FALSE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))
    print(do.call(mtext, c(paste("tSNE", title), test)))

    plotexp$logsc=TRUE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))
    print(do.call(mtext, c(paste("tSNE (Log)", title), test)))

    plotexp$logsc=FALSE; plotexp$fr = TRUE
    print(do.call(plotexpmap, c(sc, plotexp)))
    print(do.call(mtext, c(paste("FR", title), test)))

    plotexp$logsc=TRUE; plotexp$fr = TRUE
    print(do.call(plotexpmap, c(sc, plotexp)))
    print(do.call(mtext, c(paste("FR (Log)", title), test)))

    if (!is.null(plotmarkg$samples)){
        reg <- plotmarkg$samples
        plotmarkg$samples <- sub("(\\_\\d+)$","", colnames(sc@ndata))
    }
    print(do.call(plotmarkergenes, c(sc, plotmarkg)))
}

sc <- in.rdat

pdf(out.pdf)
if (perform.plotting) do.plotting(sc)
if (perform.symbolmap) do.inspect.symbolmap(sc)
if (perform.genesofinterest) do.inspect.genesofinterest(sc)
if (perform.diffgene) do.inspect.diffgene(sc)
dev.off()
