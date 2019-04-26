#!/usr/bin/Rscript --vanilla
# Density Plot Module for Galaxy
# ggplot2
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
#
# Version 1
# Cristel Thomas
#
#

library(ggplot2)
library(grid)
# Multiple plot function
# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

generateGraphFromText <- function(input, channels, output, plot_default=TRUE,
                                  flag_pdf=FALSE) {
  fcs <- read.table(input, header = TRUE, sep = "\t", check.names = FALSE)
  ## marker names
  markers <- colnames(fcs)

  if (plot_default) {
    channels <- c(grep(colnames(fcs), pattern="Forward scatter", ignore.case=TRUE),
                  grep(colnames(fcs), pattern="Side scatter", ignore.case=TRUE))
  	if (length(channels) == 0){
      channels <- c(grep(colnames(fcs), pattern="FSC"),
                    grep(colnames(fcs), pattern="SSC"))
      if (length(channels) > 2) {
        #get first FSC and corresponding SSC
        channels <- c(grep(colnames(fcs), pattern="FSC-A"),
                      grep(colnames(fcs), pattern="SSC-A"))
        if (length(channels) == 0) {
          channels <- c(grep(colnames(fcs), pattern="FSC-H"),
                        grep(colnames(fcs), pattern="SSC-H"))
          if (length(channels) == 0) {
            channels <- c(grep(colnames(fcs), pattern="FSC-W"),
                          grep(colnames(fcs), pattern="SSC-W"))
          }
        }
      }
    }
    if (length(channels) == 0) {
      warning('No forward/side scatter channels found, no plots will be generated.')
	  quit(save = "no", status = 10, runLast = FALSE)
    }
  }

  nb_markers <- length(channels)
  if (nb_markers == 1) {
    warning('There is only one marker selected to plot.')
    quit(save = "no", status = 12, runLast = FALSE)
  }
  for (j in nb_markers) {
    if (channels[j] > length(markers)){
  	  warning('Please indicate markers between 1 and ', length(markers))
  	  quit(save = "no", status = 10, runLast = FALSE)
  	}
  }

  plots <- list()
  i <- 0
  for (m in 1:(nb_markers - 1)) {
    for (n in (m+1):nb_markers) {
      x <- fcs[,channels[m]]
      y <- fcs[,channels[n]]
      df <- data.frame(x = x, y = y,
                       d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
      p <- ggplot(df) +
           geom_point(aes(x, y, col = d), size = 0.2) +
           scale_color_identity() +
           theme_bw() +
           labs(x = markers[channels[m]]) +
           labs(y = markers[channels[n]])
      i <- i + 1
      plots[[i]] <- p
    }
  }
  nb_rows <- ceiling(((nb_markers-1)*nb_markers)/4)
  h <- 400 * nb_rows
  hp <- 10 * (nb_rows/2)

  if (flag_pdf){
    pdf(output, height=hp, width=10, useDingbats=FALSE, onefile=TRUE)
      multiplot(plotlist = plots, cols = 2)
    dev.off()
  } else {
    png(output, type="cairo", width=800, height=h)
      multiplot(plotlist = plots, cols = 2)
    dev.off()
  }
}

args <- commandArgs(trailingOnly = TRUE)
channels <- list()
flag_default <- FALSE
flag_pdf <- FALSE

if (args[2]=="None" || args[2]== "" || args[2] == "i.e.:1,3,4") {
  flag_default <- TRUE
} else {
  channels <- as.numeric(strsplit(args[2], ",")[[1]])
  for (channel in channels){
    if (is.na(channel)){
      quit(save = "no", status = 11, runLast = FALSE)
    }
  }
  if (length(channels) == 1){
    warning('Please indicate more than one marker to plot.')
    quit(save = "no", status = 10, runLast = FALSE)
  }
}

if (args[4] == "PDF"){
  flag_pdf <- TRUE
}
generateGraphFromText(args[1], channels, args[3], flag_default, flag_pdf)
