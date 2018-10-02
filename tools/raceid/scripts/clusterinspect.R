#!/usr/bin/env R
VERSION = "0.1"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressPackageStartupMessages(require(RaceID))
source(args[1])

do.inspect.symbolmap <- function(sc){  
    if (!is.null(plotsym.use.typeremoveregex)){
        plotsym$types = sub(plotsym.use.typeremoveregex, "", colnames(sc@ndata))
        
        if (!is.null(plotsym.use.typeremoveregex.subselect)){
            plotsym$subset = plotsym$types[grep(plotsym.use.typeremoveregex.subselect, plotsym$types)]
        }
    }
    plotsym$fr = FALSE
    print(do.call(plotsymbolsmap, c(sc, plotsym)))
    plotsym$fr = TRUE
    print(do.call(plotsymbolsmap, c(sc, plotsym)))
}

do.inspect.genesofinterest <- function(sc){
    plotexp$logsc=FALSE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))

    plotexp$logsc=TRUE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))

    plotexp$logsc=FALSE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))

    plotexp$logsc=TRUE; plotexp$fr = FALSE
    print(do.call(plotexpmap, c(sc, plotexp)))

    if (!is.null(plotmarkg$samples)){
        reg <- plotmarkg$samples
        plotmarkg$samples <- sub("(\\_\\d+)$","", colnames(sc@ndata))
    }
    
    print(do.call(plotmarkergenes, c(sc, plotmarkg)))
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
}

sc <- in.rdat

pdf(out.pdf)
if (perform.symbolmap) do.inspect.symbolmap(sc)
if (perform.genesofinterest) do.inspect.genesofinterest(sc)
if (perform.diffgene) do.inspect.diffgene(sc)
dev.off()


