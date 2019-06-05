#!/usr/bin/env R

suppressPackageStartupMessages(require(ggplot2))

log10histoPlot <- function(title="", columncounts){
    #' Log10 histogram plot
    #'
    #' @param columncounts colSums(input_matrix)
    #' @param title Title of plot
    #' @return ggplot grob
    dfer <- data.frame(colcounts=log10(columncounts))
    p1 <- ggplot(dfer, aes(x=colcounts)) +
        geom_histogram(binwidth = 0.05, color="black",fill="white") +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title=title, y="Frequency", x="Library Size (Log10)")

    return(p1)
}

## This is calculated by the first call to contaminationPlot
## and then re-used by the second call.
ylim.max = NULL

contaminationPlot <- function(title, columndata, barcode.data, plate.data, RAW)
{
    coldata = data.frame(colsums=columndata)
    maxval = max(coldata)
    # Set once and once only
    if (is.null(ylim.max)){
        ylim.max <<- maxval + 500
    }

    drawPlates <- function(plate.data){ # 0 384 768

        plate.boundaries <- plate.data["filtered.plates",]
        if (RAW){
            plate.boundaries <- plate.data["unfilter.plates",]
        }

        plate.minyval <- -200
        plate.color <- "grey"
        plate.size <- 2
        plate.text.color <- "#7a0909"
        plate.text.alpha <- 1
        plate.text.height <- -200
        plate.height <- 4 * maxval / 5
        plate.spacing <- 35
        plate.height.text <- plate.height - 3000
        plate.text.size <- 3.5

        ## Connecting bar underneath
        gplot <<- gplot + geom_segment(aes(x=min(unname(unlist(plate.boundaries))), xend=max(unname(unlist(plate.boundaries))),
                                           y=plate.minyval, yend=plate.minyval), color=plate.color, lty=1, size=plate.size,
                                       linejoin="round", lineend="round")

        ## Overlay text
        zzz <- lapply(names(plate.boundaries), function(plate.name){
            plate.num = as.integer(sub("P","",plate.name))
            plate.xval = plate.boundaries[[plate.name]]

            ## ## If not first plate -- inner boundary, print next plate name to the right
            ## if (plate.name != names(plate.boundaries)[length(plate.boundaries)]){
            ##     gplot <<- gplot + annotate("text", x=plate.xval + plate.spacing, label=paste("Plate", plate.num + 1),
            ##                                size=plate.text.size, y=plate.height*0.7, angle=-90,
            ##                                color=plate.text.color, alpha=plate.text.alpha)
            ## }
            ## ## If not last plate -- inner boundary, print previous name to the left
            ## if (plate.name != names(plate.boundaries)[1]){
            ##     gplot <<- gplot + annotate("text", x=plate.xval - plate.spacing, label=paste("Plate", plate.num),
            ##                                size=plate.text.size, y=plate.height*0.7, angle=90,
            ##                                color=plate.text.color, alpha=plate.text.alpha)
            ## }
           
            gplot <<- gplot + geom_segment(aes(x=plate.xval, xend=plate.xval, y=plate.minyval, yend=plate.height), color=plate.color, lty=1, size=plate.size)

            if (plate.num > 0){
                ## Print under, between plates
                prev.plate.name <- paste("P", plate.num - 1, sep="")
                prev.plate.xval <- plate.boundaries[[prev.plate.name]]
                p.xval = as.integer((plate.xval + prev.plate.xval)/2)

                gplot <<- gplot + annotate("text", x=p.xval, label=paste("Plate", plate.num),
                                           size=plate.text.size, y=plate.text.height-200,
                                           color=plate.text.color, alpha=plate.text.alpha)
            }
        })

    }


    drawBatches <- function(barcode.data, plate.data){
        batch.minyval = -200
        divide.color <- "red"
        divide.style <- 2 ## dashed
        divide.size <- 0.3
        divide.text.size <- 2
        divide.height <- (maxval * 7/8) - 500
        boundary.color <- "blue"
        boundary.style <- 1
        boundary.size <- 0.5
        divide.text.true.color <- "blue"
        divide.text.false.color <- "black"
        divide.text.name.color <- "#7a0909"
        batch.text.tf.height <- divide.height * 7/8
        batch.text.tf.spacing <- 20
        batch.text.alpha <- 0.8
        batch.text.size <- 3
        batch.label.text.size <- 4
        batch.height <- maxval
        batch.text.height <- maxval * 7/8
        batch.spacing <- 24

        boundary.pos <- barcode.data["filtered_positions",]
        plate.pos <- unname(unlist(plate.data["filtered.plates",])) ## just want values

        if (RAW){
            boundary.pos <- barcode.data["unfilter_positions",]
            divider.pos <- barcode.data["filter_in_unfilter",]
            plate.pos <- unname(unlist(plate.data["unfilter.plates",])) ## just want values
        }


        zzz <- lapply(names(boundary.pos), function(batch.name){
            batch.num = as.integer(sub("B","",batch.name))
            batch.xval = unname(unlist(boundary.pos[[batch.name]]))

            # Plot boundary line only if not intersecting with a plate line
            if (!(batch.xval %in% plate.pos)){
                gplot <<- gplot + geom_segment(aes(x=batch.xval, xend=batch.xval, y=batch.minyval, yend=batch.height),
                                               color=boundary.color, lty=boundary.style, size=boundary.size)
            }

            ## Plot labels (except P0)
            if (batch.name != "B0"){
                if (RAW){
                    divide.xval = unname(unlist(divider.pos[[batch.name]]))

                    ## Plot the divider line
                    gplot <<- gplot + geom_segment(aes(x=divide.xval, xend=divide.xval, y=batch.minyval, yend=divide.height),
                                               color=divide.color, lty=divide.style, size=divide.size)

                    ## Plot the True/False labels
                    gplot <<- gplot + annotate("text", x=divide.xval - batch.text.tf.spacing, label="False positives",
                                               size=divide.text.size, y=batch.text.tf.height, angle=+90,
                                               color=divide.text.false.color, alpha=batch.text.alpha)
                    gplot <<- gplot + annotate("text", x=divide.xval + batch.text.tf.spacing, label="True positives",
                                               size=divide.text.size, y=batch.text.tf.height, angle=-90,
                                               color=divide.text.true.color, alpha=batch.text.alpha)
                    ## Plot the Batch names
                    gplot <<- gplot + annotate("text", x=divide.xval, label=batch.name,
                                               size=batch.label.text.size, y=batch.text.height, angle=0,
                                               color=divide.text.name.color, alpha=batch.text.alpha)
                } else {
                    ## Plot the Batch names between blue lines
                    ##
                    ## We calculate this by taking the middle between current batch blue lines
                    ## and previous blue line
                    prev.batch.name = paste("B", batch.num - 1, sep="")
                    prev.batch.xval = unname(unlist(boundary.pos[[prev.batch.name]]))

                    b.xval <- as.integer((batch.xval + prev.batch.xval)/2)

                    gplot <<- gplot + annotate("text", x=b.xval, label=batch.name,
                                               size=batch.label.text.size, y=batch.text.height, angle=0,
                                               color = divide.text.name.color, alpha=batch.text.alpha)
                }
            }
        })
    }


    plotCells <- function(coldata){
        gplot <<- gplot + geom_point(data=coldata, aes(x=1:nrow(coldata), y=coldata$colsums), pch = 16, cex = 1)
    }

    plotTheme <- function(){
        gplot <<- gplot + theme(plot.title = element_text(hjust = 0.5),
                                axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                axis.text.x=element_blank()) +
            labs(title=paste("Contamination Plot\n", title), y="Library Size", x="Barcode Index") +
            coord_cartesian(ylim = c(-200, ylim.max)) +
            scale_x_continuous(breaks=NULL)
    }


    gplot <- ggplot()
    plotCells(coldata)
    drawBatches(barcode.data, plate.data)
    drawPlates(plate.data)
    plotTheme()

    return(gplot)
}
