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


sanityCheck(experiment.spec, colnames(input_matrix))

                                        #
                                        #
                                        # Batch Plotting Functions
                                        #
calculatePlateIndexes <- function(plate.form, full.barcode.size){
    num.plates <<- length(names(plate.form))
    batches.per.plate <- length(plate.form[[1]])
    plate.size <- batches.per.plate * full.barcode.size
    return(seq(plate.size, num.plates * plate.size, plate.size))
}

calculateFullBarcodeIndexes <- function(num.batches, full.barcode.size){
    #' For N batches and a list of actually detected barcodes in the header,
    #' generates where the blue lines should be
    bsize <- full.barcode.size
    return(seq(bsize, num.batches * bsize, bsize))    
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


### MUST! MUST! MUST! REORDER COLUMNS SO THAT FOR EACH BATCH, IT FOLLOWS THE BARCODE ORDER GIVEN IN THE SPEC

contaminationPlot <- function(columncounts, title = "", indexes.plates, indexes.fullbc, indexes.truebc)
{
    par(mfrow=c(2,1))

    plot(columncounts, pch = 16, cex = .4, ylab = "Library Size", xlab="Cell No.",
         main = paste("Contamination Plot\n", title)
    )
    abline(v=indexes.plates, col='red', lty=3)
    abline(v=indexes.fullbc, col='blue', lty=3)
    abline(v=indexes.truebc, col='green', lty=3)

    message("Red   = plates")
    message("Blue  = across all barcodes")
    message("Green = across batch-specific barcodes ")
    message("No. cells with lib sizes > 20000 : ", sum(columncounts > 20000))

    hist(log10(columncounts), breaks = 100, xlim = c(0,5))
}

out_plot = "test.png"
png(out_plot)
contaminationPlot(
    colSums(input_matrix),
    title="pre-merge",
    calculatePlateIndexes(spec$plates, num.barcodes),
    calculateFullBarcodeIndexes(num.batches, num.barcodes),
    calculateRealBarcodeIndexes(spec$format, num.barcodes)
)
dev.off()
