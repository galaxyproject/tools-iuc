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

reorderMatrixHeaders <- function(barcodes, count.matrix, barcode.format, plate.format, sort.cells){
    #' Reorder headers to segment wanted and unwanted barcodes on opposite sides
    #' of each batch
    #'
    #' @param barcodes list of full barcodes
    #' @param count.matrix input matrix
    #' @param barcode.format batch list specifying valid barcodes for each batch
    #' @param plate.format plate list specifying plate format for each batch
    #' @param sort.cells sort cells by sizes
    #' @return list of all barcodes sorted bilaterally by batch, and true barcodes

    fixed.headers <- colnames(count.matrix)

    batch.ordering <- list()
    batch.ordering.filtered <- list()

    res <- sapply(names(barcode.format), function(key){
        rng <- as.integer(unlist(strsplit(key, '-')))
        ranges <- seq(rng[1],rng[2])

        # Barcodes wanted and unwanted for this range of batches
        barc.wanted <- barcodes[ranges]
        barc.unwant <- barcodes[!(barcodes %in% barc.wanted)]

        sub.batches <- barcode.format[[key]]  # 1,3,5,7 or 2,4,6,8

        res2 <- lapply(sub.batches, function(bat){
            batch.match <- paste("_B",bat,"_",sep="")
            headers.in.batch <- fixed.headers[grepl(batch.match, fixed.headers)]
            barcodes.in.batch <- sub(".*_([ATCGN]+)$", "\\1", headers.in.batch)

            headers.in.batch.wanted <- headers.in.batch[barcodes.in.batch %in% barc.wanted]
            headers.in.batch.unwant <- headers.in.batch[barcodes.in.batch %in% barc.unwant]

            if (sum(headers.in.batch.wanted %in% headers.in.batch.unwant) > 0){
                stop("Barcode given twice!", headers.in.batch.wanted[headers.in.batch.wanted %in% headers.in.batch.unwant])
            }

            ## Perform cell sorting if desired
            if (sort.cells){
                wanted.sorted.n <- names(sort(colSums(count.matrix[headers.in.batch.wanted])))
                unwant.sorted.n <- names(sort(colSums(count.matrix[headers.in.batch.unwant])))

                headers.in.batch.wanted <- wanted.sorted.n
                headers.in.batch.unwant <- unwant.sorted.n
            }

            ## False on the left, True on the right
            headers.in.batch.neworder <- c(headers.in.batch.unwant, headers.in.batch.wanted)
            batch.name <- paste("B", bat, sep="")

            batch.ordering[[batch.name]] <<- headers.in.batch.neworder
            batch.ordering.filtered[[batch.name]] <<- headers.in.batch.wanted
        })
    })
    ## Now we have sorted all our barcodes in each batch to the correct order
    ## we just have to sort the batches into the correct order according to plating setup
    barcode.ordering <- c()
    barcode.ordering.filtered <- c()

    res <- sapply(sortedBatchesOrPlates(names(plate.format)), function(plate.num){
        batches <- plate.format[[plate.num]]

        ## Preserve batch order on plates
        res2 <- sapply(batches, function(batch.num){
            batch.name <- paste("B", batch.num, sep="")
            barcs <- batch.ordering[[batch.name]]
            barcs.filtered <- batch.ordering.filtered[[batch.name]]

            barcode.ordering <<- c(barcode.ordering, barcs)
            barcode.ordering.filtered <<- c(barcode.ordering.filtered, barcs.filtered)
        })
    })

    return(list(
        all=barcode.ordering,
        filtered=barcode.ordering.filtered,
        debug.barcodes = batch.ordering
    ))
}