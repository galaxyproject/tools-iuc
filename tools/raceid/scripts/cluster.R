#!/usr/bin/env R
VERSION <- "0.6" # nolint

args <- commandArgs(trailingOnly = T)

if (length(args) != 1) {
     message(paste("VERSION:", VERSION))
     stop("Please provide the config file")
}

suppressWarnings(suppressPackageStartupMessages(require(RaceID)))
suppressWarnings(suppressPackageStartupMessages(require(scran)))
source(args[1])


do.filter <- function(sc) { # nolint
    if (!is.null(filt.lbatch.regexes)) {
        lar <- filt.lbatch.regexes # nolint
        nn <- colnames(sc@expdata)
        filt$LBatch <- lapply(seq_length(lar), function(m) { # nolint
            return(nn[grep(lar[[m]], nn)])
        })
    }

    sc <- do.call(filterdata, c(sc, filt)) # nolint

    ## Get histogram metrics for library size and number of features
    raw.lib <- log10(colSums(as.matrix(sc@expdata)))
    raw.feat <- log10(colSums(as.matrix(sc@expdata) > 0))
    filt.lib <- log10(colSums(getfdata(sc))) # nolint
    filt.feat <- log10(colSums(getfdata(sc) > 0)) # nolint

    if (filt.geqone) {
        filt.feat <- log10(colSums(getfdata(sc) >= 1)) # nolint
    }

    br <- 50
                                        # nolint start
    ## Determine limits on plots based on the unfiltered data
    ## (doesn't work, R rejects limits and norm data is too different
    ## to compare to exp data so let them keep their own ranges)

    ## betterrange <- function(floatval) {
    ##     return(10 * (floor(floatval / 10) + 1))
    ## }

    ## tmp.lib <- hist(raw.lib, breaks=br, plot=F)
    ## tmp.feat <- hist(raw.feat, breaks=br, plot=F)

    ## lib.y_lim <- c(0,betterrange(max(tmp.lib$counts)))
    ## lib.x_lim <- c(0,betterrange(max(tmp.lib$breaks)))

    ## feat.y_lim <- c(0,betterrange(max(tmp.feat$counts))) 
    ## feat.x_lim <- c(0,betterrange(max(tmp.feat$breaks)))
                                        # nolint end

    par(mfrow = c(2, 2))
    print(hist(raw.lib, breaks = br, main = "RawData Log10 LibSize"))
    ## , xlim=lib.x_lim, ylim=lib.y_lim)
    print(hist(raw.feat, breaks = br, main = "RawData Log10 NumFeat"))
    ##, xlim=feat.x_lim, ylim=feat.y_lim)
    print(hist(filt.lib, breaks = br, main = "FiltData Log10 LibSize"))
    ## , xlim=lib.x_lim, ylim=lib.y_lim)
    tmp <- hist(filt.feat, breaks = br, main = "FiltData Log10 NumFeat")
    ## , xlim=feat.x_lim, ylim=feat.y_lim)
    print(tmp)
    ## required, for extracting midpoint
    unq <- unique(filt.feat)
    if (length(unq) == 1) {
        abline(v = unq, col = "red", lw = 2)
        text(tmp$mids, table(filt.feat)[[1]] - 100, pos = 1,
             paste(10^unq, "\nFeatures\nin remaining\nCells", sep = ""),
             cex = 0.8)
    }

    if (filt.use.ccorrect) {
        par(mfrow = c(2, 2))
                                        # nolint start
        sc <- do.call(CCcorrect, c(sc, filt.ccc))
        print(plotdimsat(sc, change = T))
        print(plotdimsat(sc, change = F))
                                        # nolint end
    }
    return(sc)
}

do.cluster <- function(sc) { # nolint
    sc <- do.call(compdist, c(sc, clust.compdist)) # nolint
    sc <- do.call(clustexp, c(sc, clust.clustexp)) # nolint
    if (clust.clustexp$sat) {
        print(plotsaturation(sc, disp = F)) # nolint
        print(plotsaturation(sc, disp = T)) # nolint
    }
    print(plotjaccard(sc)) # nolint
    return(sc)
}

do.outlier <- function(sc) { # nolint
    sc <- do.call(findoutliers, c(sc, outlier.findoutliers)) # nolint
    if (outlier.use.randomforest) {
        sc <- do.call(rfcorrect, c(sc, outlier.rfcorrect)) # nolint
    }
    print(plotbackground(sc)) # nolint
    print(plotsensitivity(sc)) # nolint
    print(plotoutlierprobs(sc)) # nolint
    ## Heatmaps
    test1 <- list()
    test1$side <- 3
    test1$line <- 0  #1 #3

    x <- clustheatmap(sc, final = FALSE) # nolint
    print(do.call(mtext, c(paste("(Initial)"), test1)))
    x <- clustheatmap(sc, final = TRUE) # nolint
    print(do.call(mtext, c(paste("(Final)"), test1)))
    return(sc)
}

do.clustmap <- function(sc) { # nolint
    sc <- do.call(comptsne, c(sc, cluster.comptsne)) # nolint
    sc <- do.call(compfr, c(sc, cluster.compfr)) # nolint
    return(sc)
}


mkgenelist <- function(sc) {
    ## Layout
    test <- list()
    test$side <- 3
    test$line <- 0  #1 #3
    test$cex <- 0.8

    df <- c()
    options(cex = 1)
    lapply(unique(sc@cpart), function(n) {
        tryCatch({
            dg <- clustdiffgenes(sc, cl = n, pvalue = genelist.pvalue) # nolint
            dg.goi <- dg[dg$fc > genelist.foldchange, ] # nolint
            dg.goi.table <- head(dg.goi, genelist.tablelim) # nolint
            df <<- rbind(df, cbind(n, dg.goi.table))

            goi <- head(rownames(dg.goi.table), genelist.plotlim) # nolint
            print(plotmarkergenes(sc, goi)) # nolint
            buffer <- paste(rep("", 36), collapse = " ")
            print(do.call(mtext, c(paste(buffer, "Cluster ", n), test)))
            test$line <- -1
            print(do.call(mtext, c(paste(buffer, "Sig. Genes"), test)))
            test$line <- -2
            print(do.call(mtext, c(paste(buffer, "(fc > ",
                                         genelist.foldchange, ")"), test)))
        }, error = {
            print(paste("Could not print marker genes for cluster ", n))
            return(NA)
        })
    })
    write.table(df, file = out.genelist, sep = "\t", quote = F) # nolint
}


writecellassignments <- function(sc) {
    dat <- sc@cluster$kpart
    tab <- data.frame(row.names = NULL,
                      cells = names(dat),
                      cluster.initial = dat,
                      cluster.final = sc@cpart,
                      is.outlier = names(dat) %in% sc@out$out)

    write.table(tab, file = out.assignments, sep = "\t", quote = F,
                row.names = F)
}


pdf(out.pdf)

if (use.filtnormconf) {
    sc <- do.filter(sc)
    message(paste(" - Source:: genes:", nrow(sc@expdata),
                  ", cells:", ncol(sc@expdata)))
    message(paste(" - Filter:: genes:", nrow(getfdata(sc)),
                  ", cells:", ncol(getfdata(sc))))
    message(paste("         :: ",
                  sprintf("%.1f", 100 * nrow(getfdata(sc)) / nrow(sc@expdata)),
                  "% of genes remain,",
                  sprintf("%.1f", 100 * ncol(getfdata(sc)) / ncol(sc@expdata)),
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
