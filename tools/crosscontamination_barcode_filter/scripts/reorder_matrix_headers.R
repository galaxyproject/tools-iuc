#!/usr/bin/env R
##
## Reorder matrix
##
convertHeadersToSensible <- function(regex.from, regex.to, col.names){
    #' Strips headers of filenames and sets plate, batch, and barcodes
    #'
    #' @param regex.from format to extract plate, batch, and barcodes
    #' @param regex.to format to set
    #' @param matrix input matrix to rename headers
    #' @return updated names
    return(sub(regex.from, regex.to, col.names))
}

reorderMatrixHeaders <- function(barcodes, headers, barcode.format){
    #' Reorder headers to segment wanted and unwanted barcodes on opposite sides
    #' of each batch
    #'
    #' @param barcodes list of full barcodes
    #' @param headers input matrix headers
    #' @param barcode.format batch list specifying valid barcodes for each batch
    #' @return list of all barcodes sorted bilaterally by batch, and true barcodes
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
