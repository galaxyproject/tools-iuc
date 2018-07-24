#!/usr/bin/env R
## VERSION = "0.1"


## args = commandArgs(trailingOnly = T)

## if (length(args) != 1){
##     message(paste("VERSION: ", VERSION))
##     stop("Please provide the config file")
## }

## #source(args[1])

##input_matrix <- readRDS('/extra/nomansland/repos/mtekman/private/notebooks/single-cell/projects/wolfdriev_2018/source/super_matrix.rds')
input_matrix <- read.table(
    'input_matrix.tsv', stringsAsFactors = F,
    na.strings=c("NA", "-", "?", "."), header=TRUE, row.names=1
)
input_matrix[is.na(input_matrix)] <- 0

spec = list(
    barcodes = 'celseq_barcodes.192.raw',

    format = list(
        "1-96"   = c(1,3,5,7),
        "97-192" = c(2,4,6,8)
    ),
    plates = list(
        "1" = c(1,2,3,4),
        "2" = c(5,6,7,8)
    )
)

regex_to_extract_P_and_B = ".*P(\\d)_(\\d)_([ACTG]+)"
new_headers <- sub(regex_to_extract_P_and_B, "P\\1_B\\2_\\3", colnames(input_matrix))
colnames(input_matrix) <- new_headers

## At this point we expect our headers to have names: P1_B2_AACCTT
barcodes <- NULL
num.barcodes <- NULL
num.batches <- NULL
num.plates <- NULL

                                        #
                                        #
                                        # Sanity Check Methods
                                        #
checkNoMissingRanges <- function(format, barcodes){
    ranges <- c()
    res <- sapply(names(format), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        ranges <<- c(ranges, seq(rng[1],rng[2]))
    })
    full.range <- seq(min(ranges),max(ranges))
    not.in <- !(full.range %in% ranges)
    if (sum(not.in) != 0) {
        message("Warning: values[",
                paste(full.range[not.in], collapse=","),
                "] -> barcodes[",
                paste(barcodes[not.in], collapse=","),
                "] are not used."
                )
    } else {
        message("All barcodes accounted for.")
    }
    return(barcodes[!not.in])

}

assertNoMissingBarcodes <- function(headers, barcodes){
    ## headers are in format P1_B1_ACTG
    barcs.in.matrix <- unique(sort(sub(".*_.*_([ACTG]+)", "\\1", new_headers)))
    not.in <- !(barcs.in.matrix %in% barcodes)
    if (sum(not.in) > 0){
        message("Error: Barcodes in matrix not in barcodes file\n", barc.in.matrix[not.in])
        stop("")
    } else {
        message("All input matrix barcodes accounted for.")
    }
}

assertNoMissingBatches <- function(format, plates){
    batches.form = c()
    batches.plate = c()
    for (form in format){batches.form = c(batches.form, form)}
    for (plate in plates){batches.plate = c(batches.plate, plate)}

    if (length(batches.plate) != length(batches.form)){
        stop("Error: The number of batches specified in the plate do not match those given in the barcode format")
    }

    range.form <- seq(min(batches.form), max(batches.form))
    range.plate <- seq(min(batches.plate), max(batches.plate))

    if (sum(!(range.form %in% batches.form)) > 0){
        stop("Error: Missing batch in barcode format")
    }
    if (sum(!(range.plate %in% batches.plate)) > 0){
        stop("Error: Missing batch in plate format")
    }
    return(length(range.form))
}


sanityCheck <- function(spec, matrix.headers){
    barcodes <<- scan(spec$barcodes, what="", sep="\n")
    num.barcodes <<- length(barcodes)

    used.barcodes <- checkNoMissingRanges(spec$format, barcodes)
    assertNoMissingBarcodes(matrix.headers, used.barcodes)
    num.batches <<- assertNoMissingBatches(spec$format, spec$plates)
}


sanityCheck(spec, colnames(input_matrix))

                                        #
                                        #
                                        # Batch Plotting Functions
                                        #
calculatePlateIndexes <- function(plate.form, full.barcode.size){
    num.plates <<- length(names(plate.form))
    batches.per.plate <- length(plate.form[[1]])
    plate.size <- batches.per.plate * full.barcode.size
    return(seq(0, num.plates * plate.size, plate.size))
}

calculateFullBarcodeIndexes <- function(num.batches, full.barcode.size){
    #' For N batches and a list of actually detected barcodes in the header,
    #' generates where the blue lines should be
    bsize <- full.barcode.size
    return(seq(0, num.batches * bsize, bsize))
}

calculateRealBarcodeIndexes <- function(barcode.form, full.barcode.size){
    #' For N batches a list of actually USED barcodes as given by the spec,
    #' generates where the green lines should be
    batches <- c()
    res <- sapply(names(barcode.form), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        size.of.range <- length(seq(rng[1],rng[2]))
        sub.batches <- barcode.form[[key]]  # 1,3,5,7 or 2,4,6,8
        res2 <- lapply(sub.batches, function(bat){
            batches[[bat]] <<- size.of.range
        })
    })
    ## We now have sizes per batch, in order of batch
    ## Need to place these at positions after each full barcode size
    positions <- c()
    res <- sapply(1:length(batches), function(b){
        batch.start <- (b-1) * full.barcode.size
        batch.size <- batches[[b]]
        positions <<- c(positions, batch.start + batch.size)
    })
    return(positions)
}

                                        #
                                        #
                                        # Reorder matrix
                                        #
reorderMatrixHeaders <- function(headers, barcode.format){
    #' For Batch get acceptable barcodes
    form <- barcode.format
    batch.ordering <- list()
    batch.ordering.correct <- list()

    res <- sapply(names(form), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        ranges <- seq(rng[1],rng[2])

        barc.wanted <- barcodes[ranges]
        barc.unwant <- barcodes[!(barcodes %in% barc.wanted)]

        sub.batches <- form[[key]]  # 1,3,5,7 or 2,4,6,8
        res2 <- lapply(sub.batches, function(bat){
            batch_bar <- headers[grepl(paste("P\\d_B",bat,"_([ACGT]+)", sep=""), headers)]
            barcs.in.batch <- sub("P._B._([ACGT]+)", "\\1", batch_bar)
            b.wanted <- batch_bar[barcs.in.batch %in% barc.wanted]
            b.unwant <- batch_bar[barcs.in.batch %in% barc.unwant]

            if (sum(b.wanted %in% b.unwant) > 0){
                stop("Barcode given twice!", b.wanted[b.wanted %in% b.unwant])
            }
            barc_order <- c(b.wanted, b.unwant)
            batch.ordering[[bat]] <<- barc_order
            batch.ordering.correct[[bat]] <<- b.wanted
        })
    })

    barcode.ordering <- c()
    barcode.ordering.correct <- c()

    res <- lapply(1:length(batch.ordering), function(bat){
        barc_order <- batch.ordering[[bat]]
        barc_order.correct <- batch.ordering.correct[[bat]]
        barcode.ordering <<- c(barcode.ordering, barc_order)
        barcode.ordering.correct <<- c(barcode.ordering.correct, barc_order.correct)
    })

    return(list(all=barcode.ordering,correct=barcode.ordering.correct))
}

ordering <- reorderMatrixHeaders(colnames(input_matrix), spec.format)

nmatrix <- input_matrix[,ordering$all]
cmatrix <- input_matrix[,ordering$correct]

library(ggplot2)
library(gridExtra)

contaminationPlot <- function(columncounts, title = "", indexes.plates, indexes.fullbc, indexes.truebc)
{   
    par(mfrow=c(2,1))
    dfer <- data.frame(colcounts=columncounts)

    ## Remove spots where plates and full barcodes mix
    #indexes.fullbc = indexes.fullbc[!(indexes.fullbc %in% indexes.plates)]

    nit <- length(indexes.truebc)
    nif <- length(indexes.fullbc)
    nip <- length(indexes.plates)
    mval <- max(dfer)

    ## Aesthetic
    min.height <- -300
    tf.spacing.left <- 12
    tf.spacing.right <- 12
    
    tf.height <- mval - 10000
    bn.height <- mval - 2000
    plate.color <- "grey"
    plate.color.text <- "black"
    plate.alpha.text <- 0.5
    plate.height <- 2*mval/5
    
    truebcs <- data.frame(x=indexes.truebc, y=rep(min.height,nit), xend=indexes.truebc, yend=rep(mval,nit))
    fullbcs <- data.frame(x=indexes.fullbc, y=rep(min.height,nif), xend=indexes.fullbc, yend=rep(mval,nif))
    platess <- data.frame(x=indexes.plates, y=rep(min.height,nip), xend=indexes.plates, yend=rep(plate.height,nip))
    connecting.bar <- data.frame(x=min(indexes.plates), y=min.height, xend=max(indexes.plates), yend=min.height)
    

    p1 <- ggplot() +
        geom_segment(data=truebcs, aes(x=x,y=y,xend=xend,yend=yend - 4000), col='grey', lty=2, size=0.2) +
        geom_segment(data=fullbcs, aes(x=x,y=y,xend=xend,yend=yend), col='blue', lty=1, size=0.4, alpha=0.2) +
        geom_segment(data=platess, aes(x=x,y=y,xend=xend,yend=yend), col=plate.color, lty=1, size=1) +
        geom_segment(data=connecting.bar, aes(x=x,y=y + 200 ,xend=xend,yend=yend), col=plate.color, lty=1, size=1) +
        geom_point(
            date=dfer, aes(x=1:length(rownames(dfer)), y=dfer$colcounts),
            pch = 16, cex = 1) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title=paste("Contamination Plot\n", title), y="Library Size", x="Cell No.") +
        scale_y_continuous(breaks=seq(0,mval + 10000, 10000)) +
        scale_x_continuous(breaks=seq(0,length(rownames(dfer)),500))

    ## Add true/false and batch labels
    res <- lapply(indexes.truebc, function(xval){
        batch <- match(xval, indexes.truebc)
        
        p1 <<- p1 +
            annotate("text", x=xval - tf.spacing.left, size=2, y=tf.height, angle=90,
                     label=" true positives", color = "dark blue", alpha = 0.5) +
            annotate("text", x=xval + tf.spacing.right, size=2, y=tf.height, angle=-90,
                     label="false positives", color = "black", alpha = 0.5) +
            annotate("text", x=xval, size=4, y=bn.height, angle=-90,
                     label=paste("B",batch,sep=""), color = "grey", alpha = 0.8)
    })

    ## Add Plate labels
    res <- lapply(indexes.plates, function(p){
        plate.num <- match(p, indexes.plates)
        left.of <- 18
        right.of <- 18

        c.label <- paste("Plate", plate.num, sep="")
        b.label <- paste("Plate", plate.num - 1, sep="")
        siz <- 3
        hei <- plate.height - 3000

        # Right label
        if (plate.num <  length(indexes.plates)){
            p1 <<- p1 +
                annotate("text", x=p + right.of, size=siz, y=hei, angle=-90,
                         label=c.label, color = plate.color.text, alpha=plate.alpha.text)
        }

        # Left label
        if (plate.num > 1){
            p1 <<- p1 +
                annotate("text", x=p - left.of, size=siz, y=hei, angle=90,
                         label=b.label, color = plate.color.text,  alpha=plate.alpha.text)
        }
    })
        


    ggsave("test.pdf", device="pdf")

    message("Red   = plates")
    message("Blue  = across all barcodes")
    message("Green = across batch-specific barcodes ")
    message("No. cells with lib sizes > 20000 : ", sum(columncounts > 20000))

}

contaminationPlot(
    colSums(nmatrix),
    title="Pre-Filter",
    calculatePlateIndexes(spec$plates, num.barcodes),
    calculateFullBarcodeIndexes(num.batches, num.barcodes),
    calculateRealBarcodeIndexes(spec$format, num.barcodes)
)
## TODO -- filter out unwanted, and plot contamination plot and histogram as seperate pages of a PDF
##      -- also fix the axes report the minimum
