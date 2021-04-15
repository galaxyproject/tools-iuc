#!/usr/bin/env R
VERSION <- "0.5"

args <- commandArgs(trailingOnly = T)

if (length(args) != 1) {
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
## suppressWarnings(suppressPackageStartupMessages(require(scran)))  # nolint
source(args[1])


do.filter <- function(sc) {
    if (!is.null(filt.lbatch.regexes)) {
        lar <- filt.lbatch.regexes
        nn <- colnames(sc@expdata)
        filt$LBatch <- lapply(seq_len(lar), function(m) {
            return(nn[grep(lar[[m]], nn)])})
    }

    sc <- do.call(filterdata, c(sc, filt))

    ## Get histogram metrics for library size and number of features
    raw.lib <- log10(colSums(as.matrix(sc@expdata)))
    raw.feat <- log10(colSums(as.matrix(sc@expdata) > 0))
    filt.lib <- log10(colSums(as.matrix(getfdata(sc))))
    filt.feat <- log10(colSums(as.matrix(getfdata(sc) > 0)))

    if (filt.geqone) {
        filt.feat <- log10(colSums(as.matrix(getfdata(sc) >= 1)))
    }

    br <- 50
    par(mfrow = c(2, 2))
    print(hist(raw.lib, breaks = br, main = "RawData Log10 LibSize"))
    print(hist(raw.feat, breaks = br, main = "RawData Log10 NumFeat"))
    print(hist(filt.lib, breaks = br, main = "FiltData Log10 LibSize"))
    tmp <- hist(filt.feat, breaks = br, main = "FiltData Log10 NumFeat")
    print(tmp)
    ## required, for extracting midpoint
    unq <- unique(filt.feat)
    if (length(unq) == 1) {
        abline(v = unq, col = "red", lw = 2)
        text(tmp$mids, table(filt.feat)[[1]] - 100, pos = 1,
             paste(10^unq, "\nFeatures\nin remaining\nCells", sep = ""), cex = 0.8)
    }

    if (filt.use.ccorrect) {
        par(mfrow = c(2, 2))
        sc <- do.call(CCcorrect, c(sc, filt.ccc))
        print(plotdimsat(sc, change = T))
        print(plotdimsat(sc, change = F))
    }
    return(sc)
}

do.cluster <- function(sc) {
    sc <- do.call(compdist, c(sc, clust.compdist))
    sc <- do.call(clustexp, c(sc, clust.clustexp))
    if (clust.clustexp$sat) {
        print(plotsaturation(sc, disp = F))
        print(plotsaturation(sc, disp = T))
    }
    print(plotjaccard(sc))
    return(sc)
}

do.outlier <- function(sc) {
    sc <- do.call(findoutliers, c(sc, outlier.findoutliers))
    if (outlier.use.randomforest) {
        sc <- do.call(rfcorrect, c(sc, outlier.rfcorrect))
    }
    print(plotbackground(sc))
    print(plotsensitivity(sc))
    print(plotoutlierprobs(sc))
    ## Heatmaps
    test1 <- list()
    test1$side <- 3
    test1$line <- 0  #1 #3

    x <- clustheatmap(sc, final = FALSE)
    print(do.call(mtext, c(paste("(Initial)"), test1)))
    x <- clustheatmap(sc, final = TRUE)
    print(do.call(mtext, c(paste("(Final)"), test1)))
    return(sc)
}

do.clustmap <- function(sc) {
    sc <- do.call(comptsne, c(sc, cluster.comptsne))
    sc <- do.call(compfr, c(sc, cluster.compfr))
    sc <- do.call(compumap, c(sc, cluster.compumap))
    return(sc)
}


mkgenelist <- function(sc) {
    ## Layout
    test <- list()
    test$side <- 4
    test$line <- -2
    test$cex <- 0.8

    df <- c()
    options(cex = 1)
    plot.new()
    lapply(unique(sc@cpart), function(n) {
        dg <- clustdiffgenes(sc, cl = n, pvalue = genelist.pvalue)$dg

        dg.goi <- dg[dg$fc > genelist.foldchange, ]
        dg.goi.table <- head(dg.goi, genelist.tablelim)
        df <<- rbind(df, cbind(n, dg.goi.table))

        goi <- head(rownames(dg.goi.table), genelist.plotlim)

        print(plotmarkergenes(sc, goi))
        buffer <- paste(rep("", 36), collapse = " ")
        print(do.call(mtext, c(paste(buffer, "Cluster ", n), test)))
        test$line <- -1
        print(do.call(mtext, c(paste(buffer, "Sig. Genes"), test)))
        test$line <- 0
        print(do.call(mtext, c(paste(buffer, "(fc > ",
                                     genelist.foldchange, ")"), test)))
    })
    write.table(df, file = out.genelist, sep = "\t", quote = F)
}


writecellassignments <- function(sc) {
    dat <- sc@cluster$kpart
    tab <- data.frame(row.names = NULL,
                      cells = names(dat),
                      cluster.initial = dat,
                      cluster.final = sc@cpart,
                      is.outlier = names(dat) %in% sc@out$out)

    write.table(tab, file = out.assignments, sep = "\t",
                quote = F, row.names = F)
}


pdf(out.pdf)

if (use.filtnormconf) {
    sc <- do.filter(sc)
    message(paste(" - Source:: genes:", nrow(sc@expdata),
                  ", cells:", ncol(sc@expdata)))
    message(paste(" - Filter:: genes:", nrow(as.matrix(getfdata(sc))),
                  ", cells:", ncol(as.matrix(getfdata(sc)))))
    message(paste("         :: ",
                  sprintf("%.1f", 100 * nrow(as.matrix(
                                            getfdata(sc))) / nrow(sc@expdata)),
                  "% of genes remain, ",
                  sprintf("%.1f", 100 * ncol(as.matrix(
                                            getfdata(sc))) / ncol(sc@expdata)),
                  "% of cells remain"))
    write.table(as.matrix(sc@ndata), file = out.table, col.names = NA,
                row.names = T, sep = "\t", quote = F)
}

if (use.cluster) {
    par(mfrow = c(2, 2))
    sc <- do.cluster(sc)

    par(mfrow = c(2, 2))
    sc <- do.outlier(sc)

    par(mfrow = c(2, 2), mar = c(1, 1, 6, 1))
    sc <- do.clustmap(sc)

    mkgenelist(sc)
    writecellassignments(sc)
}

dev.off()

saveRDS(sc, out.rdat)
