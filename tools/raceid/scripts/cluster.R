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

    if (filt.use.ccorrect){
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
    x <- clustheatmap(sc, final=FALSE)
    x <- clustheatmap(sc, final=TRUE)
    return(sc)
}

do.clustmap <- function(sc){
    sc <- do.call(comptsne, c(sc, cluster.comptsne))
    sc <- do.call(compfr, c(sc, cluster.compfr))
    return(sc)
}


mkgenelist <- function(sc){
    df <- c()
    lapply(unique(sc@cpart), function(n){
        dg <- clustdiffgenes(sc, cl=n, pvalue=genelist.pvalue)

        dg.goi <- dg[dg$fc > genelist.foldchange,]
        dg.goi.table <- head(dg.goi, genelist.tablelim)
        df <<- rbind(df, cbind(n, dg.goi.table))

        goi <- head(rownames(dg.goi.table), genelist.plotlim)
        print(plotmarkergenes(sc, goi))
    })
    write.table(df, file=out.genelist, sep="\t", quote=F)
}

pdf(out.pdf)
par(mfrow=c(2,2))

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
