#!/usr/bin/Rscript --vanilla
######################################################################
#                  Copyright (c) 2016 Northrop Grumman.
#                          All rights reserved.
######################################################################
# ImmPort FCS conversion program
# Authors: Yue Liu and Yu "Max" Qian
#
# Reference: FCSTrans: An open source software system for FCS
#            file conversion and data transformation
#            Qian Y, Liu Y, Campbell J, Thomson E, Kong YM,
#            Scheuermann RH. 2012 Cytometry Part A. 81A(5)
#            doi.org/10.1002/cyto.a.22037
#
# To run in R
# 1) library(flowCore)
# 2) source("FCSTrans.R")
# 3) transformFCS("filename")
#
#
# Automated Gating of Lymphocytes with FlowDensity
# Authors of FlowDensity: Jafar Taghiyar, Mehrnoush Malek
#
# Reference: flowDensity: reproducing manual gating of flow
#            cytometry data by automated density-based cell
#            population identification
#            Malek M, Taghiyar MJ, Chong L, Finak G,
#            Gottardo R, Brinkman RR. 2015 Bioinformatics 31(4)
#            doi: 10.1093/bioinformatics/btu677
#
#
# Version 1.5
# March 2016 -- added lines to run directly from command line (cristel thomas)
# May 2016 -- added automated gating (cristel thomas)
# August 2016 -- added options for data transformation (cristel thomas)
# April 2017 -- added logicle to transformation options (cristel thomas)
# July 2017 -- added options for outputs (cristel thomas)

library(flowCore)
library(flowDensity)
library(GEOmap)
#
# Set output to 0 when input is less than cutoff value
#
ipfloor <- function (x, cutoff=0, target=0) {
  y <- x
  if (x <= cutoff) {
    y <- target
  }
  return(y)
}
#
# Set output to 0 when input is less than cutoff value
#
ipceil <- function (x, cutoff=0, target=0) {
  y <- x
  if (x >= cutoff) {
    y <- target
  }
  return(y)
}
#
# Calculation core of iplogicle
#
iplogicore <- function (x, w, r, d, scale) {
  tol <- .Machine$double.eps^0.8
  maxit <- as.integer(5000)
  d <- d * log(10)
  scale <- scale / d
  p <- if (w == 0) {
         1
       } else {
         uniroot(function(p) -w + 2 * p * log(p)/(p + 1), c(.Machine$double.eps,
         2 * (w + d)))$root
       }
  a <- r * exp(-(d - w))
  b <- 1
  c <- r * exp(-(d - w)) * p^2
  d <- 1/p
  f <- a * (p^2 - 1)
  y <- .Call("flowCore_biexponential_transform", PACKAGE= 'flowCore',
             as.double(x), a, b, c, d, f, w, tol, maxit)
  y <- sapply(y * scale, ipfloor)
  return(y)
}
#
# Function for calculating w
#
iplogiclew <- function (w, cutoff=-111, r=262144, d=4.5, scale=1) {
  if (w > d)
    w <- d
  y <- iplogicore(cutoff, w, r, d, scale) - .Machine$double.eps^0.6
  return(y)
}
#
# ImmPort logicle function - convert fluorescent marker values to channel output
#
iplogicle <- function (x, r=262144, d=4.5, range=4096, cutoff=-111, w=-1) {
  if (w > d) {
    stop("Negative range decades must be smaller than total number of decades")
  }
  if (w < 0) {
    w = uniroot(iplogiclew, c(0, d), cutoff=cutoff)$root
  }
  y <- iplogicore(x, w, r, d, range)
  return(y)
}
#
# Convert fluorescent values to channel output using log transformation
#
iplog <- function(x) {
  x <- sapply(x, ipfloor, cutoff=1, target=1)
  y <- 1024 * log10(x) - 488.6
  return(y)
}
#
# ImmPort linear function - convert scatter values to channel output
# linear transformation
#
ipscatter <- function (x, channelrange=262144) {
  y <- 4095.0 * x / channelrange
  y <- sapply(y, ipfloor)
  y <- sapply(y, ipceil, cutoff=4095, target=4095)
  return(y)
}
#
# ImmPort time function - convert time values to channel output
# linear transformation
iptime <- function (x, channelrange) {
  # use simple cutoff for now
  y <- sapply(x, ipfloor)
  return(y)
}
#
# Determine the type of marker. Marker type is used
# to determine type of transformation to apply for this channel.
# Before 2010 FLUO_AREA type used iplogicile and
# FLOU_NON_AREA type used iplog. In 2010 Yue, changed code so
# all fluorescent channels use iplogicle. Below is the note from SVN
#
# Version 1.1
# 2010-07-02
# -----------
# Added data type checking on both FCS version 2 and 3
# Removed log conversion for non-area fluorescent channel
# Applied logicle conversion for all fluorescent channels
#
# The GenePattern version uses iplog for FLOU_NON_AREA, rather
# than iplogicle.
#
getMarkerType <- function(name,debug=FALSE) {
  type <- ""
  prefix2 <- toupper(substr(name, 1, 2))
  prefix3 <- toupper(substr(name, 1, 3))
  prefix4 <- toupper(substr(name, 1, 4))
  if (prefix2 == "FS" || prefix2 == "SS") {
    type <- "SCATTER"
  } else if (prefix3 == "FSC" || prefix3 == "SSC") {
    type <- "SCATTER"
  } else if (prefix4 == "TIME") {
    type <- "TIME"
  } else {
    pieces <- unlist(strsplit(name, "-"))
    if (toupper(pieces[length(pieces)]) == "A") {
      type <- "FLUO_AREA"
    } else {
      type <- "FLUO_NON_AREA"
    }
  }
  if (debug) {
    print(paste("Marker:", name, ", Type:", type))
  }
  return(type)
}
#
# Scale data
#
scaleData <- function(data, channelrange=0) {
  datamax <- range(data)[2]  # range() returns [min, max]
  if (datamax > channelrange) {
    channelrange <- datamax
  }
  #if (channelrange == 0) {
  #    channelrange = range(data)[2]
  #}
  data <- 262144 * data / channelrange
  return(data)
}
#
# Check if AccuriData. Accuri data needs different conversion
#
isAccuriData <- function(keywords) {
  isTRUE(as.character(keywords$"$CYT") == "Accuri C6")
}
#
# Convert FCS file
#
convertFCS <- function(fcs, debug=FALSE) {
  # Check file type and FCS version
  if (class(fcs)[1] != "flowFrame") {
    print("convertFCS requires flowFrame object as input")
    return(FALSE)
  }
  keywords <- keyword(fcs)
  markers <- colnames(fcs)
  params <- fcs@parameters
  list_description <- fcs@description

  if (debug) {
    print("****Inside convertFCS")
    print(paste("    FCS version:", keywords$FCSversion))
    print(paste("    DATATYPE:", keywords['$DATATYPE']))
  }
  if (keywords$FCSversion == "2" || keywords$FCSversion == "3" ||
      keywords$FCSversion == "3.1" ) {
    datatype <- unlist(keywords['$DATATYPE'])
    if (datatype == 'F') {
      # Apply compensation if available and requested
      # Process fcs expression data, using transformation
      # based on category of the marker.
      fcs_exprs <- exprs(fcs)
      fcs_channel <- NULL
      for (i in 1:length(markers)){
        markertype <- getMarkerType(markers[i], debug)
        rangekeyword <- paste("$P", i, "R", sep="")
        flowcore_min <- paste("flowCore_", rangekeyword, "min", sep="")
        flowcore_max <- paste("flowCore_", rangekeyword, "max", sep="")
        channelrange <- as.numeric(keywords[rangekeyword])
        if (debug) {
          print(paste("    Marker name:", markers[i]))
          print(paste("    Marker type:", markertype))
          print(paste("    Range value:", keywords[rangekeyword]))
        }

        if (markertype == "TIME") {
          channel <- iptime(fcs_exprs[, i])
        } else {
          if (markertype == "SCATTER") {
            channel <- ipscatter(scaleData(fcs_exprs[, i], channelrange))
          } else {
            # Apply logicle transformation on fluorescent channels
            channel <- iplogicle(scaleData(fcs_exprs[, i], channelrange))
          }
          # adjust range in parameters and list description
          if (params@data$range[i] > 4096){
            params@data$range[i] <- 4096
            params@data$minRange[i] <- 0
            params@data$maxRange[i] <- 4096
            list_description[rangekeyword] <- 4096
            list_description[flowcore_min] <- 0
            list_description[flowcore_max] <- 4096
          }
        }
        fcs_channel <- cbind(fcs_channel, round(channel))
      }
      colnames(fcs_channel) <- markers
    } else {
      if (datatype != 'I') {
        print(paste("Data type", datatype, "in FCS 3 is not supported"))
      }
      fcs_channel <- exprs(fcs)
      colnames(fcs_channel) <- markers
    }
  } else {
    print(paste("FCS version", keyword(fcs)$FCSversion, "is not supported"))
    fcs_channel <- exprs(fcs)
    colnames(fcs_channel) <- markers
  }
  newfcs <- flowFrame(fcs_channel, params, list_description)
  return(newfcs)
}
#
# Starting function for processing a FCS file
#
processFCSFile <- function(input_file, output_file="", compensate=FALSE,
                           outformat="flowtext", gate=FALSE,
                           graph_file="", report="", method="",
                           scaling_factor, logicle_w=0.5, logicle_t=262144,
                           logicle_m=4.5, debug=FALSE) {
  #
  # Generate the file names for the output_file
  #
  pieces <- unlist(strsplit(input_file, .Platform$file.sep))
  filename <- pieces[length(pieces)]
  if (debug) {
    print (paste("Converting file: ",input_file))
    print (paste("Original file name: ",filename))
    print (paste("Output file name: ",output_file))
  }
  fcs <- read.FCS(input_file, transformation=F)
  keywords <- keyword(fcs)
  markers <- colnames(fcs)
  print_markers <- as.vector(pData(parameters(fcs))$desc)
  # Update print_markers if the $P?S not in the FCS file
  for (i in 1:length(print_markers)) {
    if (is.na(print_markers[i])) {
      print_markers[[i]] <- markers[i]
    }
  }
  #
  # Compensate
  #
  spill <- keywords$SPILL

  if (is.null(spill) == FALSE && compensate == TRUE) {
    if (debug) {
      print("Attempting compensation")
    }
    tryCatch({fcs = compensate(fcs, spill)},
              error = function(ex) {str(ex); })
  }
  #
  # Transform the data
  #
  transformed_data <- fcs
  channels_to_exclude <- c(grep(colnames(fcs), pattern="FSC"),
                           grep(colnames(fcs), pattern="SSC"),
                           grep(colnames(fcs), pattern="Time"))
  list_channels <- colnames(fcs)[-channels_to_exclude]
  if (isAccuriData(keywords)) {
    print("Accuri data is not supported")
  } else if (method != "None"){
    if (method == "fcstrans"){
      transformed_data <- convertFCS(fcs, debug)
    } else if (method == "logicle_auto"){
      lgcl <- estimateLogicle(fcs, channels = list_channels)
      transformed_data <- transform(fcs, lgcl)
    } else {
      if (method == "arcsinh"){
        trans <- arcsinhTransform(transformationId="defaultArcsinhTransform",
                                  a = 0, b = scaling_factor, c = 0)
      } else  if (method == "logicle"){
        trans <- logicleTransform(w = logicle_w, t = logicle_t, m = logicle_m)
      }
      translist <- transformList(list_channels, trans)
      transformed_data <- transform(fcs, translist)
    }
  }
  trans_gated_data <- transformed_data
  #
  # Gate data
  #
  if (gate){
    # check that there are SSC and FSC channels to gate on
    chans <- c(grep(colnames(transformed_data), pattern="FSC"),
               grep(colnames(transformed_data), pattern="SSC"))
    totalchans <- chans
    if (length(chans) > 2) {
      #get first FSC and corresponding SSC
      chans <- c(grep(colnames(transformed_data), pattern="FSC-A"),
                 grep(colnames(transformed_data), pattern="SSC-A"))
      if (length(chans) == 0) {
        chans <- c(grep(colnames(transformed_data), pattern="FSC-H"),
                   grep(colnames(transformed_data), pattern="SSC-H"))
        if (length(chans) == 0) {
          chans <- c(grep(colnames(transformed_data), pattern="FSC-W"),
                     grep(colnames(transformed_data), pattern="SSC-W"))
        }
      }
    }
    if (length(chans) == 0) {
      warning('No forward/side scatter channels found, gating aborted.')
    } else {
      # gate lymphocytes
      lymph <- flowDensity(obj=transformed_data, channels=chans,
                           position=c(TRUE, NA),
                           debris.gate=c(TRUE, FALSE))
      # gate singlets if A and H/W
      if (length(totalchans) > 2) {
          trans_gated_data <- getflowFrame(flowDensity(obj=lymph,
                                           singlet.gate=TRUE))
      } else {
        trans_gated_data <- getflowFrame(lymph)
      }
      # report
      pregating_summary <- capture.output(summary(transformed_data))
      pregating_dim <- capture.output(dim(transformed_data))
      postgating_summary <- capture.output(summary(trans_gated_data))
      postgating_dim <- capture.output(dim(trans_gated_data))
      sink(report)
      cat("#########################\n")
      cat("##    BEFORE GATING    ##\n")
      cat("#########################\n")
      cat(pregating_dim, pregating_summary, sep="\n")
      cat("\n#########################\n")
      cat("##    AFTER  GATING    ##\n")
      cat("#########################\n")
      cat(postgating_dim, postgating_summary, sep="\n")
      sink()
      # plots
      time_channel <- grep(toupper(colnames(transformed_data)), pattern="TIME")
      nb_markers <- length(colnames(transformed_data)) - length(time_channel)
      nb_rows <- ceiling(((nb_markers-1)*nb_markers)/4)
      h <- 400 * nb_rows
      maxrange <- transformed_data@parameters@data$range[1]

      png(graph_file, type="cairo", height=h, width=800)
      par(mfrow=c(nb_rows,2))
      for (m in 1:(nb_markers - 1)) {
        for (n in (m+1):nb_markers) {
          plotDens(transformed_data, c(m,n), xlab = print_markers[m],
                   ylab = print_markers[n], main = "Before Gating",
                   ylim = c(0, maxrange), xlim = c(0, maxrange))
          plotDens(trans_gated_data, c(m,n), xlab = print_markers[m],
                   ylab = print_markers[n], main = "After Gating",
                   ylim = c(0, maxrange), xlim = c(0, maxrange))
        }
      }
      dev.off()
    }
  }
  if (outformat=="FCS") {
    write.FCS(trans_gated_data, output_file)
  } else if (outformat=="flowFrame") {
    saveRDS(trans_gated_data, file = output_file)
  } else {
    output_data <- exprs(trans_gated_data)
    colnames(output_data) <- print_markers
    write.table(output_data, file=output_file, quote=F,
                row.names=F,col.names=T, sep='\t', append=F)
  }
}
# Convert FCS file using FCSTrans logicile transformation
# @param input_file     FCS file to be transformed
# @param output_file    FCS file transformed ".txt" extension
# @param compensate     Flag indicating whether to apply compensation
#                       matrix if it exists.
transformFCS <- function(input_file, output_file, compensate=FALSE,
                         outformat="flowtext", gate=FALSE, graph_file="",
                         report_file="", trans_met="", scaling_factor=1/150,
                         w=0.5, t=262144, m=4.5, debug=FALSE) {
  isValid <- F
  # Check file beginning matches FCS standard
  tryCatch({
    isValid <- isFCSfile(input_file)
  }, error = function(ex) {
    print (paste("    ! Error in isFCSfile", ex))
  })
  if (isValid) {
    processFCSFile(input_file, output_file, compensate, outformat,
                   gate, graph_file, report_file, trans_met, scaling_factor,
                   w, t, m)
  } else {
    print (paste(input_file, "does not meet FCS standard"))
  }
}
#
# Run FCS Gate-Trans
#
args <- commandArgs(trailingOnly = TRUE)
graphs <- ""
report <- ""
gate <- FALSE
trans_method <- "None"
scaling_factor <- 1 / 150
w <- 0.5
t <- 262144
m <- 4.5
if (args[5]!="None") {
  gate <- TRUE
  graphs <- args[5]
  report <- args[6]
}
if (args[7]!="None"){
  trans_method <- args[7]
  if (args[7] == "arcsinh"){
    scaling_factor <- 1 / as.numeric(args[8])
  } else if (args[7] == "logicle"){
    w <- args[8]
    t <- args[9]
    m <- args[10]
  }
}
transformFCS(args[1], args[2], args[3], args[4], gate, graphs,
             report, trans_method, scaling_factor, w, t, m)
