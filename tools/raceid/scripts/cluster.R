#!/usr/bin/env R
VERSION = "0.2"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
suppressWarnings(suppressPackageStartupMessages(require(scran)))
source(args[1])


do.filter <- function(sc){
    if (!is.null(filt.lbatch.regexes)){
        lar <- filt.lbatch.regexes
        nn <- colnames(sc@expdata)
        filt$LBatch <- lapply(1:length(lar), function(m){ return( nn[grep(lar[[m]], nn)] ) })
    }

    sc <- do.call(filterdata, c(sc, filt))

    ## Get histogram metrics for library size and number of features
    raw.lib <- log10(colSums(as.matrix(sc@expdata)))
    raw.feat <- log10(rowSums(as.matrix(sc@expdata)>0))
    filt.lib <- log10(colSums(getfdata(sc)))
    filt.feat <- log10(rowSums(getfdata(sc)>0))

    br <- 50
    ## Determine limits on plots based on the unfiltered data
    ## (doesn't work, R rejects limits and norm data is too different to compare to exp data
    ##  so let them keep their own ranges)

    ## betterrange <- function(floatval){
    ##     return(10 * (floor(floatval / 10) + 1))
    ## }

    ## tmp.lib <- hist(raw.lib, breaks=br, plot=F)
    ## tmp.feat <- hist(raw.feat, breaks=br, plot=F)

    ## lib.y_lim <- c(0,betterrange(max(tmp.lib$counts)))
    ## lib.x_lim <- c(0,betterrange(max(tmp.lib$breaks)))

    ## feat.y_lim <- c(0,betterrange(max(tmp.feat$counts)))
    ## feat.x_lim <- c(0,betterrange(max(tmp.feat$breaks)))

    par(mfrow=c(2,2))
    print(hist(raw.lib, breaks=br, main="RawData Log10(LibSize)")) # , xlim=lib.x_lim, ylim=lib.y_lim)
    print(hist(raw.feat, breaks=br, main="RawData Log10(NumFeat)")) #, xlim=feat.x_lim, ylim=feat.y_lim)
    print(hist(filt.lib, breaks=br, main="FiltData Log10(LibSize)")) # , xlim=lib.x_lim, ylim=lib.y_lim)
    tmp <- hist(filt.feat, breaks=br, main="FiltData Log10(NumFeat)") # , xlim=feat.x_lim, ylim=feat.y_lim)
    print(tmp) # required, for extracting midpoint
    unq <- unique(filt.feat)
    if (length(unq) == 1){
        text(tmp$mids, table(filt.feat)[[1]] - 100, pos=1, paste(format(10^unq, scientific=T, digits=3),
                                                           " Features in all Cells", sep=""), cex=0.8)
    }


    if (filt.use.ccorrect){
        par(mfrow=c(2,2))
        sc <- do.call(CCcorrect, c(sc, filt.ccc))
        print(plotdimsat(sc, change=T))
        print(plotdimsat(sc, change=F))
    }
    return(sc)
}

do.cluster <- function(sc){
    sc <- do.call(compdist, c(sc, clust.compdist))
    sc <- do.call(clustexp, c(sc, clust.clustexp))
    if (clust.clustexp$sat){
        print(plotsaturation(sc, disp=F))
        print(plotsaturation(sc, disp=T))
    }
    print(plotjaccard(sc))
    return(sc)
}

do.outlier <- function(sc){
    sc <- do.call(findoutliers, c(sc, outlier.findoutliers))
    if (outlier.use.randomforest){
        sc <- do.call(rfcorrect, c(sc, outlier.rfcorrect))
    }
    print(plotbackground(sc))
    print(plotsensitivity(sc))
    print(plotoutlierprobs(sc))
    ## Heatmaps
    test1 <- list()
    test1$side = 3
    test1$line = 0  #1 #3

    x <- clustheatmap(sc, final=FALSE)
    print(do.call(mtext, c(paste("(Initial)"), test1)))  ## spacing is a hack
    x <- clustheatmap(sc, final=TRUE)
    print(do.call(mtext, c(paste("(Final)"), test1)))  ## spacing is a hack
    return(sc)
}

do.clustmap <- function(sc){
    sc <- do.call(comptsne, c(sc, cluster.comptsne))
    sc <- do.call(compfr, c(sc, cluster.compfr))
    return(sc)
}


mkgenelist <- function(sc){
    ## Layout
    test <- list()
    test$side = 3
    test$line = 0  #1 #3
    test$cex = 0.8

    df <- c()
    options(cex = 1)
    lapply(unique(sc@cpart), function(n){
        dg <- clustdiffgenes(sc, cl=n, pvalue=genelist.pvalue)

        dg.goi <- dg[dg$fc > genelist.foldchange,]
        dg.goi.table <- head(dg.goi, genelist.tablelim)
        df <<- rbind(df, cbind(n, dg.goi.table))

        goi <- head(rownames(dg.goi.table), genelist.plotlim)
        print(plotmarkergenes(sc, goi))
        print(do.call(mtext, c(paste("                               Cluster ",n), test)))  ## spacing is a hack
        test$line=-1
        print(do.call(mtext, c(paste("                               Sig. Genes"), test)))  ## spacing is a hack
        test$line=-2
        print(do.call(mtext, c(paste("                               (fc > ", genelist.foldchange,")"), test)))  ## spacing is a hack

    })
    write.table(df, file=out.genelist, sep="\t", quote=F)
}

pdf(out.pdf)

if (use.filtnormconf){
    sc <- do.filter(sc)
    message(paste(" - Source:: genes:",nrow(sc@expdata),", cells:",ncol(sc@expdata)))
    message(paste(" - Filter:: genes:",nrow(sc@ndata),", cells:",ncol(sc@ndata)))
    message(paste("         :: ",
                  sprintf("%.1f", 100 * nrow(sc@ndata)/nrow(sc@expdata)), "% of genes remain,",
                  sprintf("%.1f", 100 * ncol(sc@ndata)/ncol(sc@expdata)), "% of cells remain"))
}

if (use.cluster){
    par(mfrow=c(2,2))
    sc <- do.cluster(sc)

    par(mfrow=c(2,2))
    sc <- do.outlier(sc)

    par(mfrow=c(2,2), mar=c(1,1,6,1))
    sc <- do.clustmap(sc)

    mkgenelist(sc)
}

dev.off()

saveRDS(sc, out.rdat)
