#!/usr/bin/env R

require(ggplot2)

log10histoPlot <- function(columncounts, title=""){
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

contaminationPlot <- function(columncounts, title = "",
                              indexes.plates, indexes.fullbc, indexes.truebc,
                              filtered=FALSE)
{
    #' Plots true and false barcodes
    #'
    #'
    #' @param columncounts colSums(input_matrix)
    #' @param title plot title
    #' @param indexes.plates plate line positions
    #' @param indexes.fullbc full batch line positions
    #' @param indexes.truebc true batch line positions
    #' @param filtered specifies whether the positions have been adjusted for true barcodes
    #' @return ggplot grob
    dfer <- data.frame(colcounts=columncounts)

    ## Remove indexes where plates and full barcodes mix
    indexes.fullbc = indexes.fullbc[!(indexes.fullbc %in% indexes.plates)]
    
    nit <- length(indexes.truebc)
    nif <- length(indexes.fullbc)
    nip <- length(indexes.plates)
    mval <- max(dfer)
  
    ## Aesthetics
    min.height <- -300
    tf.spacing.left <- 12
    tf.spacing.right <- 12   
    tf.height <- mval - 10000
    bn.height <- mval - 2000
    plate.color <- "grey"
    plate.text.color <- "black"
    plate.text.alpha <- 0.5
    plate.text.size <- 3
    plate.height <- 2* mval / 5
    plate.spacing <- if (filtered) 12 else 24
    plate.height.text <- plate.height - 3000
    
    truebcs <- data.frame(x=indexes.truebc, y=rep(min.height,nit), xend=indexes.truebc, yend=rep(mval,nit))
    fullbcs <- data.frame(x=indexes.fullbc, y=rep(min.height,nif), xend=indexes.fullbc, yend=rep(mval,nif))
    platess <- data.frame(x=indexes.plates, y=rep(min.height,nip), xend=indexes.plates, yend=rep(plate.height,nip))
    connecting.bar <- data.frame(x=min(indexes.plates), y=min.height, xend=max(indexes.plates), yend=min.height)  
      
    p1 <- ggplot()

    if (!filtered){
        p1 <- p1 +
            geom_segment(data=truebcs, aes(x=x,y=y,xend=xend,yend=yend - 4000), col='grey', lty=2, size=0.2) + 
            geom_segment(data=fullbcs, aes(x=x,y=y,xend=xend,yend=yend), col='blue', lty=1, size=0.4, alpha=0.2)
    }
    else {
        p1 <- p1 +
            geom_segment(data=truebcs, aes(x=x,y=y,xend=xend,yend=yend), col='blue', lty=1, size=0.4, alpha=0.2)
    }
    
    p1 <- p1 +
        geom_segment(data=platess, aes(x=x,y=y,xend=xend,yend=yend), col=plate.color, lty=1, size=1) +
        geom_segment(data=connecting.bar, aes(x=x,y=y + 200 ,xend=xend,yend=yend), col=plate.color, lty=1, size=1) +
        geom_point(
            data=dfer, aes(x=1:length(rownames(dfer)), y=dfer$colcounts),
            pch = 16, cex = 1) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
              axis.text.x=element_blank()) +
        labs(title=paste("Contamination Plot\n", title), y="Library Size", x="Barcode Index") +
        scale_y_continuous(breaks=seq(0,mval + 10000, 10000)) +
        scale_x_continuous(breaks=NULL)

    ## Add true/false and batch labels
    res <- lapply(indexes.truebc, function(xval){
        batch <- match(xval, indexes.truebc)

        if (!filtered){
            p1 <<- p1 +
                annotate("text", x=xval - tf.spacing.left, size=2, y=tf.height, angle=90, label=" true positives", color = "dark blue", alpha = 0.5) +
                annotate("text", x=xval + tf.spacing.right, size=2, y=tf.height, angle=-90,label="false positives", color = "black", alpha = 0.5) +
                annotate("text", x=xval , size=4, y=bn.height, angle=-90, label=paste("B",batch,sep=""), color = "grey", alpha = 0.8)
        }
        else {
            p1 <<- p1 +
                annotate("text", x=xval - 48, size=4, y=bn.height, angle=-90, label=paste("B",batch,sep=""), color = "grey", alpha = 0.8)
        }
    })

    ## Add Plate labels
    res <- lapply(indexes.plates, function(p){
        plate.num <- match(p, indexes.plates)
        c.label <- paste("Plate", plate.num, sep="")
        b.label <- paste("Plate", plate.num - 1, sep="")

        # Right label
        if (plate.num <  length(indexes.plates)){
            p1 <<- p1 +
                annotate("text", x=p + plate.spacing, size=plate.text.size, y=plate.height.text, angle=-90,
                         label=c.label, color = plate.text.color, alpha=plate.text.alpha)
        }

        # Left label
        if (plate.num > 1){
            p1 <<- p1 +
                annotate("text", x=p - plate.spacing, size=plate.text.size, y=plate.height.text, angle=90,
                         label=b.label, color = plate.text.color,  alpha=plate.text.alpha)
        }
    })
        
    return(p1)
}


