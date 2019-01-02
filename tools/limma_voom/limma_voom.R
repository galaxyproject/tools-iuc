# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: htmlPath", "R", 1, "character"      -Path to html file linking to other outputs
#       outPath", "o", 1, "character"       -Path to folder to write all output to
#       filesPath", "j", 2, "character"     -JSON list object if multiple files input
#       matrixPath", "m", 2, "character"    -Path to count matrix
#       factFile", "f", 2, "character"      -Path to factor information file
#       factInput", "i", 2, "character"     -String containing factors if manually input
#       annoPath", "a", 2, "character"      -Path to input containing gene annotations
#       contrastData", "C", 1, "character"  -String containing contrasts of interest
#       cpmReq", "c", 2, "double"           -Float specifying cpm requirement
#       cntReq", "z", 2, "integer"          -Integer specifying minimum total count requirement
#       sampleReq", "s", 2, "integer"       -Integer specifying cpm requirement
#       normCounts", "x", 0, "logical"      -String specifying if normalised counts should be output
#       rdaOpt", "r", 0, "logical"          -String specifying if RData should be output
#       lfcReq", "l", 1, "double"           -Float specifying the log-fold-change requirement
#       pValReq", "p", 1, "double"          -Float specifying the p-value requirement
#       pAdjOpt", "d", 1, "character"       -String specifying the p-value adjustment method
#       normOpt", "n", 1, "character"       -String specifying type of normalisation used
#       robOpt", "b", 0, "logical"          -String specifying if robust options should be used
#       trend", "t", 1, "double"            -Float for prior.count if limma-trend is used instead of voom
#       weightOpt", "w", 0, "logical"       -String specifying if voomWithQualityWeights should be used
#       topgenes", "G", 1, "integer"        -Integer specifying no. of genes to highlight in volcano and heatmap
#       treatOpt", "T", 0, "logical"        -String specifying if TREAT function should be used
#       plots, "P", 1, "character"          -String specifying additional plots to be created
#
# OUT:
#       Density Plots (if filtering)
#       Box Plots (if normalising)
#       MDS Plot
#       Voom/SA plot
#       MD Plot
#       Volcano Plot
#       Heatmap
#       Expression Table
#       HTML file linking to the ouputs
# Optional:
#       Normalised counts Table
#       RData file
#
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Jun 2017, Jan 2018, May 2018

# Record starting time
timeStart <- as.character(Sys.time())

# Load all required libraries
library(methods, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)
library(scales, quietly=TRUE, warn.conflicts=FALSE)
library(getopt, quietly=TRUE, warn.conflicts=FALSE)
library(gplots, quietly=TRUE, warn.conflicts=FALSE)

################################################################################
### Function Declaration
################################################################################
# Function to sanitise contrast equations so there are no whitespaces
# surrounding the arithmetic operators, leading or trailing whitespace
sanitiseEquation <- function(equation) {
    equation <- gsub(" *[+] *", "+", equation)
    equation <- gsub(" *[-] *", "-", equation)
    equation <- gsub(" *[/] *", "/", equation)
    equation <- gsub(" *[*] *", "*", equation)
    equation <- gsub("^\\s+|\\s+$", "", equation)
    return(equation)
}

# Function to sanitise group information
sanitiseGroups <- function(string) {
    string <- gsub(" *[,] *", ",", string)
    string <- gsub("^\\s+|\\s+$", "", string)
    return(string)
}

# Function to change periods to whitespace in a string
unmake.names <- function(string) {
    string <- gsub(".", " ", string, fixed=TRUE)
    return(string)
}

# Generate output folder and paths
makeOut <- function(filename) {
    return(paste0(opt$outPath, "/", filename))
}

# Generating design information
pasteListName <- function(string) {
    return(paste0("factors$", string))
}

# Create cata function: default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)
cata <- function(..., file = opt$htmlPath, sep = "", fill = FALSE, labels = NULL,
               append = TRUE) {
    if (is.character(file))
        if (file == "")
            file <- stdout()
        else if (substring(file, 1L, 1L) == "|") {
            file <- pipe(substring(file, 2L), "w")
            on.exit(close(file))
        }
        else {
        file <- file(file, ifelse(append, "a", "w"))
      on.exit(close(file))
        }
    .Internal(cat(list(...), file, sep, fill, labels, append))
}

# Function to write code for html head and title
HtmlHead <- function(title) {
    cata("<head>\n")
    cata("<title>", title, "</title>\n")
    cata("</head>\n")
}

# Function to write code for html links
HtmlLink <- function(address, label=address) {
    cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

# Function to write code for html images
HtmlImage <- function(source, label=source, height=500, width=500) {
    cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
    cata("\" width=\"", width, "\"/>\n")
}

# Function to write code for html list items
ListItem <- function(...) {
    cata("<li>", ..., "</li>\n")
}

TableItem <- function(...) {
    cata("<td>", ..., "</td>\n")
}

TableHeadItem <- function(...) {
    cata("<th>", ..., "</th>\n")
}

################################################################################
### Input Processing
################################################################################

# Collect arguments from command line
args <- commandArgs(trailingOnly=TRUE)

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    "htmlPath", "R", 1, "character",
    "outPath", "o", 1, "character",
    "filesPath", "j", 2, "character",
    "matrixPath", "m", 2, "character",
    "factFile", "f", 2, "character",
    "factInput", "i", 2, "character",
    "annoPath", "a", 2, "character",
    "contrastData", "C", 1, "character",
    "cpmReq", "c", 1, "double",
    "totReq", "y", 0, "logical",
    "cntReq", "z", 1, "integer",
    "sampleReq", "s", 1, "integer",
    "filtCounts", "F", 0, "logical",
    "normCounts", "x", 0, "logical",
    "rdaOpt", "r", 0, "logical",
    "lfcReq", "l", 1, "double",
    "pValReq", "p", 1, "double",
    "pAdjOpt", "d", 1, "character",
    "normOpt", "n", 1, "character",
    "robOpt", "b", 0, "logical",
    "trend", "t", 1, "double",
    "weightOpt", "w", 0, "logical",
    "topgenes", "G", 1, "integer",
    "treatOpt", "T", 0, "logical",
    "plots", "P", 1, "character",
    "libinfoOpt", "L", 0, "logical"),
    byrow=TRUE, ncol=4)
opt <- getopt(spec)


if (is.null(opt$matrixPath) & is.null(opt$filesPath)) {
    cat("A counts matrix (or a set of counts files) is required.\n")
    q(status=1)
}

if (is.null(opt$cpmReq)) {
    filtCPM <- FALSE
} else {
    filtCPM <- TRUE
}

if (is.null(opt$cntReq) || is.null(opt$sampleReq)) {
    filtSmpCount <- FALSE
} else {
    filtSmpCount <- TRUE
}

if (is.null(opt$totReq)) {
    filtTotCount <- FALSE
} else {
    filtTotCount <- TRUE
}

if (is.null(opt$rdaOpt)) {
    wantRda <- FALSE
} else {
    wantRda <- TRUE
}

if (is.null(opt$annoPath)) {
    haveAnno <- FALSE
} else {
    haveAnno <- TRUE
}

if (is.null(opt$filtCounts)) {
    wantFilt <- FALSE
} else {
    wantFilt <- TRUE
}

if (is.null(opt$normCounts)) {
    wantNorm <- FALSE
} else {
    wantNorm <- TRUE
}

if (is.null(opt$robOpt)) {
    wantRobust <- FALSE
} else {
    wantRobust <- TRUE
}

if (is.null(opt$weightOpt)) {
    wantWeight <- FALSE
} else {
    wantWeight <- TRUE
}

if (is.null(opt$trend)) {
    wantTrend <- FALSE
    deMethod <- "limma-voom"
} else {
    wantTrend <- TRUE
    deMethod <- "limma-trend"
    priorCount <- opt$trend
}

if (is.null(opt$treatOpt)) {
    wantTreat <- FALSE
} else {
    wantTreat <- TRUE
}

if (is.null(opt$libinfoOpt)) {
    wantLibinfo <- FALSE
} else {
    wantLibinfo <- TRUE
}


if (!is.null(opt$filesPath)) {
    # Process the separate count files (adapted from DESeq2 wrapper)
    library("rjson")
    parser <- newJSONParser()
    parser$addData(opt$filesPath)
    factorList <- parser$getObject()
    factors <- sapply(factorList, function(x) x[[1]])
    filenamesIn <- unname(unlist(factorList[[1]][[2]]))
    sampleTable <- data.frame(sample=basename(filenamesIn),
                            filename=filenamesIn,
                            row.names=filenamesIn,
                            stringsAsFactors=FALSE)
    for (factor in factorList) {
        factorName <- factor[[1]]
        sampleTable[[factorName]] <- character(nrow(sampleTable))
        lvls <- sapply(factor[[2]], function(x) names(x))
        for (i in seq_along(factor[[2]])) {
            files <- factor[[2]][[i]][[1]]
            sampleTable[files,factorName] <- lvls[i]
        }
        sampleTable[[factorName]] <- factor(sampleTable[[factorName]], levels=lvls)
    }
    rownames(sampleTable) <- sampleTable$sample
    rem <- c("sample","filename")
    factors <- sampleTable[, !(names(sampleTable) %in% rem), drop=FALSE]

    #read in count files and create single table
    countfiles <- lapply(sampleTable$filename, function(x){read.delim(x, row.names=1)})
    counts <- do.call("cbind", countfiles)

} else {
    # Process the single count matrix
    counts <- read.table(opt$matrixPath, header=TRUE, sep="\t", strip.white=TRUE, stringsAsFactors=FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[ , -1]
    countsRows <- nrow(counts)

    # Process factors
    if (is.null(opt$factInput)) {
            factorData <- read.table(opt$factFile, header=TRUE, sep="\t", strip.white=TRUE)
            # order samples as in counts matrix
            factorData <- factorData[match(colnames(counts), factorData[, 1]), ]
            factors <- factorData[, -1, drop=FALSE]
    }  else {
            factors <- unlist(strsplit(opt$factInput, "|", fixed=TRUE))
            factorData <- list()
            for (fact in factors) {
                newFact <- unlist(strsplit(fact, split="::"))
                factorData <- rbind(factorData, newFact)
            } # Factors have the form: FACT_NAME::LEVEL,LEVEL,LEVEL,LEVEL,... The first factor is the Primary Factor.

            # Set the row names to be the name of the factor and delete first row
            row.names(factorData) <- factorData[, 1]
            factorData <- factorData[, -1]
            factorData <- sapply(factorData, sanitiseGroups)
            factorData <- sapply(factorData, strsplit, split=",")
            factorData <- sapply(factorData, make.names)
            # Transform factor data into data frame of R factor objects
            factors <- data.frame(factorData)
    }
}

 # if annotation file provided
if (haveAnno) {
    geneanno <- read.table(opt$annoPath, header=TRUE, sep="\t", quote= "", strip.white=TRUE, stringsAsFactors=FALSE)
}

#Create output directory
dir.create(opt$outPath, showWarnings=FALSE)

# Split up contrasts seperated by comma into a vector then sanitise
contrastData <- unlist(strsplit(opt$contrastData, split=","))
contrastData <- sanitiseEquation(contrastData)
contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)

plots <- character()
if (!is.null(opt$plots)) {
    plots <- unlist(strsplit(opt$plots, split=","))
}

denOutPng <- makeOut("densityplots.png")
denOutPdf <- makeOut("densityplots.pdf")
cpmOutPdf <- makeOut("cpmplots.pdf")
boxOutPng <- makeOut("boxplots.png")
boxOutPdf <- makeOut("boxplots.pdf")
mdsscreeOutPng <- makeOut("mdsscree.png")
mdsscreeOutPdf <- makeOut("mdsscree.pdf")
mdsxOutPdf <- makeOut("mdsplot_extra.pdf")
mdsxOutPng <- makeOut("mdsplot_extra.png")
mdsamOutPdf <- makeOut("mdplots_samples.pdf")
mdOutPdf <- character() # Initialise character vector
volOutPdf <- character()
heatOutPdf <- character()
stripOutPdf <- character()
mdvolOutPng <- character()
topOut <- character()
glimmaOut <- character()
for (i in 1:length(contrastData)) {
    con <- contrastData[i]
    con <- gsub("\\(|\\)", "", con)
    mdOutPdf[i] <- makeOut(paste0("mdplot_", con, ".pdf"))
    volOutPdf[i] <- makeOut(paste0("volplot_", con, ".pdf"))
    heatOutPdf[i] <- makeOut(paste0("heatmap_", con, ".pdf"))
    stripOutPdf[i] <- makeOut(paste0("stripcharts_", con, ".pdf"))
    mdvolOutPng[i] <- makeOut(paste0("mdvolplot_", con, ".png"))
    topOut[i] <- makeOut(paste0(deMethod, "_", con, ".tsv"))
    glimmaOut[i] <- makeOut(paste0("glimma_", con, "/MD-Plot.html"))
}
filtOut <- makeOut(paste0(deMethod, "_", "filtcounts"))
normOut <- makeOut(paste0(deMethod, "_", "normcounts"))
rdaOut <- makeOut(paste0(deMethod, "_analysis.RData"))
sessionOut <- makeOut("session_info.txt")

# Initialise data for html links and images, data frame with columns Label and
# Link
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
                        stringsAsFactors=FALSE)

# Initialise vectors for storage of up/down/neutral regulated counts
upCount <- numeric()
downCount <- numeric()
flatCount <- numeric()

################################################################################
### Data Processing
################################################################################

# Extract counts and annotation data
print("Extracting counts")
data <- list()
data$counts <- counts
if (haveAnno) {
  # order annotation by genes in counts (assumes gene ids are in 1st column of geneanno)
  annoord <- geneanno[match(row.names(counts), geneanno[,1]), ]
  data$genes <- annoord
} else {
  data$genes <- data.frame(GeneID=row.names(counts))
}

# Creating naming data
samplenames <- colnames(data$counts)
sampleanno <- data.frame("sampleID"=samplenames, factors)

# Creating colours for the groups
cols <- as.numeric(factors[, 1])
col.group <- palette()[cols]

# If filter crieteria set, filter out genes that do not have a required cpm/counts in a required number of
# samples. Default is no filtering
preFilterCount <- nrow(data$counts)
nsamples <- ncol(data$counts)

if (filtCPM || filtSmpCount || filtTotCount) {

    if (filtTotCount) {
        keep <- rowSums(data$counts) >= opt$cntReq
    } else if (filtSmpCount) {
        keep <- rowSums(data$counts >= opt$cntReq) >= opt$sampleReq
    } else if (filtCPM) {
        myCPM <- cpm(data$counts)
        thresh <- myCPM >= opt$cpmReq 
        keep <- rowSums(thresh) >= opt$sampleReq

        if ("c" %in% plots) {
            # Plot CPM vs raw counts (to check threshold)
            pdf(cpmOutPdf, width=6.5, height=10)
            par(mfrow=c(3, 2))
            for (i in 1:nsamples) {
                plot(data$counts[, i], myCPM[, i], xlim=c(0,50), ylim=c(0,3), main=samplenames[i], xlab="Raw counts", ylab="CPM")
                abline(v=10, col="red", lty=2, lwd=2)
                abline(h=opt$cpmReq, col=4)
            }
            linkName <- "CpmPlots.pdf"
            linkAddr <- "cpmplots.pdf"
            linkData <- rbind(linkData, data.frame(Label=linkName, Link=linkAddr, stringsAsFactors=FALSE))
            invisible(dev.off())
        }
    }

    data$counts <- data$counts[keep, ]
    data$genes <- data$genes[keep, , drop=FALSE]

    if (wantFilt) {
        print("Outputting filtered counts")
        filt_counts <- data.frame(data$genes, data$counts)
        write.table(filt_counts, file=filtOut, row.names=FALSE, sep="\t", quote=FALSE)
        linkData <- rbind(linkData, data.frame(Label=paste0(deMethod, "_", "filtcounts.tsv"), Link=paste0(deMethod, "_", "filtcounts"), stringsAsFactors=FALSE))
    }

    # Plot Density
    if ("d" %in% plots) {
        # PNG
        png(denOutPng, width=1000, height=500)
        par(mfrow=c(1,2), cex.axis=0.8)

        # before filtering
        lcpm1 <- cpm(counts, log=TRUE)
        plot(density(lcpm1[, 1]), col=col.group[1], lwd=2, las=2, main="", xlab="")
        title(main="Density Plot: Raw counts", xlab="Log-cpm")
        for (i in 2:nsamples){
            den <- density(lcpm1[, i])
            lines(den$x, den$y, col=col.group[i], lwd=2)
        }

        # after filtering
        lcpm2 <- cpm(data$counts, log=TRUE)
        plot(density(lcpm2[,1]), col=col.group[1], lwd=2, las=2, main="", xlab="")
        title(main="Density Plot: Filtered counts", xlab="Log-cpm")
        for (i in 2:nsamples){
            den <- density(lcpm2[, i])
            lines(den$x, den$y, col=col.group[i], lwd=2)
        }
        legend("topright", samplenames, text.col=col.group, bty="n")
        imgName <- "Densityplots.png"
        imgAddr <- "densityplots.png"
        imageData <- rbind(imageData, data.frame(Label=imgName, Link=imgAddr, stringsAsFactors=FALSE))
        invisible(dev.off())

        # PDF
        pdf(denOutPdf, width=14)
        par(mfrow=c(1,2), cex.axis=0.8)
        plot(density(lcpm1[, 1]), col=col.group[1], lwd=2, las=2, main="", xlab="")
        title(main="Density Plot: Raw counts", xlab="Log-cpm")
        for (i in 2:nsamples){
            den <- density(lcpm1[, i])
            lines(den$x, den$y, col=col.group[i], lwd=2)
        }
        plot(density(lcpm2[, 1]), col=col.group[1], lwd=2, las=2, main="", xlab="")
        title(main="Density Plot: Filtered counts", xlab="Log-cpm")
        for (i in 2:nsamples){
            den <- density(lcpm2[, i])
            lines(den$x, den$y, col=col.group[i], lwd=2)
        }
        legend("topright", samplenames, text.col=col.group, bty="n")
        linkName <- "DensityPlots.pdf"
        linkAddr <- "densityplots.pdf"
        linkData <- rbind(linkData, data.frame(Label=linkName, Link=linkAddr, stringsAsFactors=FALSE))
        invisible(dev.off())
    }
}

postFilterCount <- nrow(data$counts)
filteredCount <- preFilterCount-postFilterCount

# Generating the DGEList object "y"
print("Generating DGEList object")
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
y <- new("DGEList", data)

print("Generating Design")
# Name rows of factors according to their sample
row.names(factors) <- names(data$counts)
factorList <- sapply(names(factors), pasteListName)
formula <- "~0"
for (i in 1:length(factorList)) {
    formula <- paste(formula,factorList[i], sep="+")
}
formula <- formula(formula)
design <- model.matrix(formula)
for (i in 1:length(factorList)) {
    colnames(design) <- gsub(factorList[i], "", colnames(design), fixed=TRUE)
}

# Calculating normalising factors
print("Calculating Normalisation Factors")
logcounts <- y #store for plots
y <- calcNormFactors(y, method=opt$normOpt)

# Generate contrasts information
print("Generating Contrasts")
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

################################################################################
### Data Output
################################################################################

# Plot Box plots (before and after normalisation)
if (opt$normOpt != "none" & "b" %in% plots) {
    png(boxOutPng, width=1000, height=500)
    par(mfrow=c(1,2), mar=c(6,4,2,2)+0.1)
    labels <- colnames(counts)

    lcpm1 <- cpm(y$counts, log=TRUE)
    boxplot(lcpm1, las=2, col=col.group, xaxt="n", xlab="")
    axis(1, at=seq_along(labels), labels = FALSE)
    abline(h=median(lcpm1), col=4)
    text(x=seq_along(labels), y=par("usr")[3]-1, srt=45, adj=1, labels=labels, xpd=TRUE)
    title(main="Box Plot: Unnormalised counts", ylab="Log-cpm")

    lcpm2 <- cpm(y, log=TRUE)
    boxplot(lcpm2, las=2, col=col.group, xaxt="n",  xlab="")
    axis(1, at=seq_along(labels), labels = FALSE)
    text(x=seq_along(labels), y=par("usr")[3]-1, srt=45, adj=1, labels=labels, xpd=TRUE)
    abline(h=median(lcpm2), col=4)
    title(main="Box Plot: Normalised counts", ylab="Log-cpm")

    imgName <- "Boxplots.png"
    imgAddr <- "boxplots.png"
    imageData <- rbind(imageData, data.frame(Label=imgName, Link=imgAddr, stringsAsFactors=FALSE))
    invisible(dev.off())

    pdf(boxOutPdf, width=14)
    par(mfrow=c(1,2), mar=c(6,4,2,2)+0.1)
    boxplot(lcpm1, las=2, col=col.group, xaxt="n", xlab="")
    axis(1, at=seq_along(labels), labels = FALSE)
    abline(h=median(lcpm1), col=4)
    text(x=seq_along(labels), y=par("usr")[3]-1, srt=45, adj=1, labels=labels, xpd=TRUE)
    title(main="Box Plot: Unnormalised counts", ylab="Log-cpm")
    boxplot(lcpm2, las=2, col=col.group, xaxt="n",  xlab="")
    axis(1, at=seq_along(labels), labels = FALSE)
    text(x=seq_along(labels), y=par("usr")[3]-1, srt=45, adj=1, labels=labels, xpd=TRUE)
    abline(h=median(lcpm2), col=4)
    title(main="Box Plot: Normalised counts", ylab="Log-cpm")
    linkName <- "BoxPlots.pdf"
    linkAddr <- "boxplots.pdf"
    linkData <- rbind(linkData, data.frame(Label=linkName, Link=linkAddr, stringsAsFactors=FALSE))
    invisible(dev.off())
}

# Plot MDS
print("Generating MDS plot")
labels <- names(counts)

# Scree plot (Variance Explained) code copied from Glimma

# get column of matrix
getCols <- function(x, inds) {
  x[, inds, drop=FALSE]
}

x <- cpm(y, log=TRUE)
ndim <- nsamples - 1
nprobes <- nrow(x)
top <- 500
top <- min(top, nprobes)
cn <- colnames(x)
bad <- rowSums(is.finite(x)) < nsamples

if (any(bad)) {
  warning("Rows containing infinite values have been removed")
  x <- x[!bad, , drop=FALSE]
}

dd <- matrix(0, nrow=nsamples, ncol=nsamples, dimnames=list(cn, cn))
topindex <- nprobes - top + 1L
for (i in 2L:(nsamples)) {
  for (j in 1L:(i - 1L)) {
    dists <- (getCols(x, i) - getCols(x, j))^2
    dists <- sort.int(dists, partial = topindex )
    topdist <- dists[topindex:nprobes]
    dd[i, j] <- sqrt(mean(topdist))
  }
}

a1 <- suppressWarnings(cmdscale(as.dist(dd), k=min(ndim, 8), eig=TRUE))
eigen <- data.frame(name = 1:min(ndim, 8), eigen = round(a1$eig[1:min(ndim, 8)]/sum(a1$eig), 2))

png(mdsscreeOutPng, width=1000, height=500)
par(mfrow=c(1, 2))
plotMDS(y, labels=samplenames, col=as.numeric(factors[, 1]), main="MDS Plot: Dims 1 and 2")
barplot(eigen$eigen, names.arg=eigen$name,  main = "Scree Plot: Variance Explained", xlab = "Dimension", ylab = "Proportion", las=1)
imgName <- paste0("MDSPlot_", names(factors)[1], ".png")
imgAddr <- "mdsscree.png"
imageData <- rbind(imageData, data.frame(Label=imgName, Link=imgAddr, stringsAsFactors=FALSE))
invisible(dev.off())

pdf(mdsscreeOutPdf, width=14)
par(mfrow=c(1, 2))
plotMDS(y, labels=samplenames, col=as.numeric(factors[, 1]), main="MDS Plot: Dims 1 and 2")
barplot(eigen$eigen, names.arg=eigen$name,  main = "Scree Plot: Variance Explained", xlab = "Dimension", ylab = "Proportion", las=1)
linkName <- paste0("MDSPlot_", names(factors)[1], ".pdf")
linkAddr <- "mdsscree.pdf"
linkData <- rbind(linkData, data.frame(Label=linkName, Link=linkAddr, stringsAsFactors=FALSE))
invisible(dev.off())

if ("x" %in% plots) {
    png(mdsxOutPng, width=1000, height=500)
    par(mfrow=c(1, 2))
    for (i in 2:3) {
        dim1 <- i
        dim2 <- i + 1
        plotMDS(y, dim=c(dim1, dim2), labels=samplenames, col=as.numeric(factors[, 1]), main=paste("MDS Plot: Dims", dim1, "and", dim2))
    }
    imgName <- paste0("MDSPlot_extra.png")
    imgAddr <- paste0("mdsplot_extra.png")
    imageData <- rbind(imageData, data.frame(Label=imgName, Link=imgAddr, stringsAsFactors=FALSE))
    invisible(dev.off())

    pdf(mdsxOutPdf, width=14)
    par(mfrow=c(1, 2))
    for (i in 2:3) {
        dim1 <- i
        dim2 <- i + 1
        plotMDS(y, dim=c(dim1, dim2), labels=samplenames, col=as.numeric(factors[, 1]), main=paste("MDS Plot: Dims", dim1, "and", dim2))
    }
    linkName <- "MDSPlot_extra.pdf"
    linkAddr <- "mdsplot_extra.pdf"
    linkData <- rbind(linkData, data.frame(Label=linkName, Link=linkAddr, stringsAsFactors=FALSE))
    invisible(dev.off())
}

if ("m" %in% plots) {
    # Plot MD plots for individual samples
    print("Generating MD plots for samples")
    pdf(mdsamOutPdf, width=6.5, height=10)
    par(mfrow=c(3, 2))
    for (i in 1:nsamples) {
        if (opt$normOpt != "none") {
            plotMD(logcounts, column=i, main=paste(colnames(logcounts)[i], "(before)"))
            abline(h=0, col="red", lty=2, lwd=2)
        }
        plotMD(y, column=i)
        abline(h=0, col="red", lty=2, lwd=2)
    }
    linkName <- "MDPlots_Samples.pdf"
    linkAddr <- "mdplots_samples.pdf"
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())
}


if (wantTrend) {
    # limma-trend approach
    logCPM <- cpm(y, log=TRUE, prior.count=opt$trend)
    fit <- lmFit(logCPM, design)
    fit$genes <- y$genes
    fit <- contrasts.fit(fit, contrasts)
    if (wantRobust) {
        fit <- eBayes(fit, trend=TRUE, robust=TRUE)
    } else {
        fit <- eBayes(fit, trend=TRUE, robust=FALSE)
    }
    # plot fit with plotSA
    saOutPng <- makeOut("saplot.png")
    saOutPdf <- makeOut("saplot.pdf")

    png(saOutPng, width=500, height=500)
    plotSA(fit, main="SA Plot")
    imgName <- "SAPlot.png"
    imgAddr <- "saplot.png"
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())

    pdf(saOutPdf, width=14)
    plotSA(fit, main="SA Plot")
    linkName <- "SAPlot.pdf"
    linkAddr <- "saplot.pdf"
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())

    plotData <- logCPM

    # Save normalised counts (log2cpm)
    if (wantNorm) {
        write.table(logCPM, file=normOut, row.names=TRUE, sep="\t", quote=FALSE)
        linkData <- rbind(linkData, c((paste0(deMethod, "_", "normcounts.tsv")), (paste0(deMethod, "_", "normcounts"))))
    }
} else {
    # limma-voom approach
    voomOutPdf <- makeOut("voomplot.pdf")
    voomOutPng <- makeOut("voomplot.png")

    if (wantWeight) {
        # Creating voom data object and plot
        png(voomOutPng, width=1000, height=500)
        vData <- voomWithQualityWeights(y, design=design, plot=TRUE)
        imgName <- "VoomPlot.png"
        imgAddr <- "voomplot.png"
        imageData <- rbind(imageData, c(imgName, imgAddr))
        invisible(dev.off())

        pdf(voomOutPdf, width=14)
        vData <- voomWithQualityWeights(y, design=design, plot=TRUE)
        linkName <- "VoomPlot.pdf"
        linkAddr <- "voomplot.pdf"
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())

        # Generating fit data and top table with weights
        wts <- vData$weights
        voomFit <- lmFit(vData, design, weights=wts)

    } else {
        # Creating voom data object and plot
        png(voomOutPng, width=500, height=500)
        vData <- voom(y, design=design, plot=TRUE)
        imgName <- "VoomPlot"
        imgAddr <- "voomplot.png"
        imageData <- rbind(imageData, c(imgName, imgAddr))
        invisible(dev.off())

        pdf(voomOutPdf)
        vData <- voom(y, design=design, plot=TRUE)
        linkName <- "VoomPlot.pdf"
        linkAddr <- "voomplot.pdf"
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())

        # Generate voom fit
        voomFit <- lmFit(vData, design)
    }

     # Save normalised counts (log2cpm)
    if (wantNorm) {
        norm_counts <- data.frame(vData$genes, vData$E)
        write.table(norm_counts, file=normOut, row.names=FALSE, sep="\t", quote=FALSE)
        linkData <- rbind(linkData, c((paste0(deMethod, "_", "normcounts.tsv")), (paste0(deMethod, "_", "normcounts"))))
    }

    # Fit linear model and estimate dispersion with eBayes
    voomFit <- contrasts.fit(voomFit, contrasts)
    if (wantRobust) {
        fit <- eBayes(voomFit, robust=TRUE)
    } else {
        fit <- eBayes(voomFit, robust=FALSE)
    }
    plotData <- vData
}

 # Save library size info
if (wantLibinfo) {
    efflibsize <- round(y$samples$lib.size * y$samples$norm.factors)
    libsizeinfo <- cbind(y$samples, EffectiveLibrarySize=efflibsize)
    libsizeinfo$lib.size <- round(libsizeinfo$lib.size)
    names(libsizeinfo)[names(libsizeinfo)=="sampleID"] <- "SampleID"
    names(libsizeinfo)[names(libsizeinfo)=="lib.size"] <- "LibrarySize"
    names(libsizeinfo)[names(libsizeinfo)=="norm.factors"] <- "NormalisationFactor"
    write.table(libsizeinfo, file="libsizeinfo", row.names=FALSE, sep="\t", quote=FALSE)
}

print("Generating DE results")

if (wantTreat) {
    print("Applying TREAT method")
    if (wantRobust) {
        fit <- treat(fit, lfc=opt$lfcReq, robust=TRUE)
    } else {
        fit <- treat(fit, lfc=opt$lfcReq, robust=FALSE)
    }
}

status = decideTests(fit, adjust.method=opt$pAdjOpt, p.value=opt$pValReq,
                       lfc=opt$lfcReq)
sumStatus <- summary(status)

for (i in 1:length(contrastData)) {
    con <- contrastData[i]
    con <- gsub("\\(|\\)", "", con)
    # Collect counts for differential expression
    upCount[i] <- sumStatus["Up", i]
    downCount[i] <- sumStatus["Down", i]
    flatCount[i] <- sumStatus["NotSig", i]

    # Write top expressions table
    if (wantTreat) {
        top <- topTreat(fit, coef=i, number=Inf, sort.by="P")
    } else{
        top <- topTable(fit, coef=i, number=Inf, sort.by="P")
    }
    write.table(top, file=topOut[i], row.names=FALSE, sep="\t", quote=FALSE)
    linkName <- paste0(deMethod, "_", con, ".tsv")
    linkAddr <- paste0(deMethod, "_", con, ".tsv")
    linkData <- rbind(linkData, c(linkName, linkAddr))

    # Plot MD (log ratios vs mean average) using limma package on weighted
    pdf(mdOutPdf[i])
    limma::plotMD(fit, status=status[, i], coef=i,
        main=paste("MD Plot:", unmake.names(con)),
        hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1),
        xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)
    linkName <- paste0("MDPlot_", con, ".pdf")
    linkAddr <- paste0("mdplot_", con, ".pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())

    # Generate Glimma interactive MD plot and table, requires annotation file (assumes gene labels/symbols in 2nd column)
    if (haveAnno) {
        # make gene labels unique to handle NAs
        geneanno <- y$genes
        geneanno[, 2] <- make.unique(geneanno[, 2])
        Glimma::glMDPlot(fit, coef=i, counts=y$counts, anno=geneanno, groups=factors[, 1],
             status=status[, i], sample.cols=col.group,
             main=paste("MD Plot:", unmake.names(con)), side.main=colnames(y$genes)[2],
             folder=paste0("glimma_", unmake.names(con)), launch=FALSE)
        linkName <- paste0("Glimma_MDPlot_", con, ".html")
        linkAddr <- paste0("glimma_", con, "/MD-Plot.html")
        linkData <- rbind(linkData, c(linkName, linkAddr))
    }

    # Plot Volcano
    pdf(volOutPdf[i])
    if (haveAnno) {
        # labels must be in second column currently
        labels <- fit$genes[, 2]
    } else {
        labels <- fit$genes$GeneID
    }
    limma::volcanoplot(fit, coef=i,
        main=paste("Volcano Plot:", unmake.names(con)),
        highlight=opt$topgenes,
        names=labels)
    linkName <- paste0("VolcanoPlot_", con, ".pdf")
    linkAddr <- paste0("volplot_", con, ".pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())

    # PNG of MD and Volcano
    png(mdvolOutPng[i], width=1000, height=500)
    par(mfrow=c(1, 2), mar=c(5,4,2,2)+0.1, oma=c(0,0,3,0))

    # MD plot
    limma::plotMD(fit, status=status[, i], coef=i, main="MD Plot",
        hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1),
        xlab="Average Expression", ylab="logFC")
    abline(h=0, col="grey", lty=2)

    # Volcano
    if (haveAnno) {
        # labels must be in second column currently
        limma::volcanoplot(fit, coef=i, main="Volcano Plot",
            highlight=opt$topgenes,
            names=fit$genes[, 2])
    } else {
        limma::volcanoplot(fit, coef=i, main="Volcano Plot",
            highlight=opt$topgenes,
            names=fit$genes$GeneID)
    }

    imgName <- paste0("MDVolPlot_", con)
    imgAddr <- paste0("mdvolplot_", con, ".png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    title(paste0("Contrast: ", unmake.names(con)), outer=TRUE, cex.main=1.5)
    invisible(dev.off())

    if ("h" %in% plots) {
        # Plot Heatmap
        topgenes <- rownames(top[1:opt$topgenes, ])
        if (wantTrend) {
            topexp <- plotData[topgenes, ]
        } else {
            topexp <- plotData$E[topgenes, ]
        }
        pdf(heatOutPdf[i])
        mycol <- colorpanel(1000,"blue","white","red")
        if (haveAnno) {
            # labels must be in second column currently
            labels <- top[topgenes, 2]
        } else {
            labels <- rownames(topexp)
        }
        heatmap.2(topexp, scale="row", Colv=FALSE, Rowv=FALSE, dendrogram="none",
            main=paste("Contrast:", unmake.names(con), "\nTop", opt$topgenes, "genes by adj.P.Val"),
            trace="none", density.info="none", lhei=c(2,10), margin=c(8, 6), labRow=labels, cexRow=0.7, srtCol=45,
            col=mycol, ColSideColors=col.group)
        linkName <- paste0("Heatmap_", con, ".pdf")
        linkAddr <- paste0("heatmap_", con, ".pdf")
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())
    }

    if ("s" %in% plots) {
        # Plot Stripcharts of top genes
        pdf(stripOutPdf[i], title=paste("Contrast:", unmake.names(con)))
        par(mfrow = c(3,2), cex.main=0.8, cex.axis=0.8)
        cols <- unique(col.group)

        for (j in 1:length(topgenes)) {
            lfc <- round(top[topgenes[j], "logFC"], 2)
            pval <- round(top[topgenes[j], "adj.P.Val"], 5)
            if (wantTrend) {
                stripchart(plotData[topgenes[j], ] ~ factors[, 1], vertical=TRUE, las=2, pch=16, cex=0.8, cex.lab=0.8, col=cols,
                    method="jitter", ylab="Normalised log2 expression", main=paste0(labels[j], "\nlogFC=", lfc, ", adj.P.Val=", pval))
            } else {
                stripchart(plotData$E[topgenes[j], ] ~ factors[, 1], vertical=TRUE, las=2, pch=16, cex=0.8, cex.lab=0.8, col=cols, 
                    method="jitter", ylab="Normalised log2 expression", main=paste0(labels[j], "\nlogFC=", lfc, ", adj.P.Val=", pval))
            }
        }
        linkName <- paste0("Stripcharts_", con, ".pdf")
        linkAddr <- paste0("stripcharts_", con, ".pdf")
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())
    }
}
sigDiff <- data.frame(Up=upCount, Flat=flatCount, Down=downCount)
row.names(sigDiff) <- contrastData

# Save relevant items as rda object
if (wantRda) {
    print("Saving RData")
    if (wantWeight) {
      save(counts, data, y, status, plotData, labels, factors, wts, fit, top, contrastData, contrasts, design,
           file=rdaOut, ascii=TRUE)
    } else {
      save(counts, data, y, status, plotData, labels, factors, fit, top, contrastData, contrasts, design,
           file=rdaOut, ascii=TRUE)
    }
    linkData <- rbind(linkData, c((paste0(deMethod, "_analysis.RData")), (paste0(deMethod, "_analysis.RData"))))
}

# Record session info
writeLines(capture.output(sessionInfo()), sessionOut)
linkData <- rbind(linkData, c("Session Info", "session_info.txt"))

# Record ending time and calculate total run time
timeEnd <- as.character(Sys.time())
timeTaken <- capture.output(round(difftime(timeEnd,timeStart), digits=3))
timeTaken <- gsub("Time difference of ", "", timeTaken, fixed=TRUE)
################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file=opt$htmlPath)

cata("<html>\n")

cata("<body>\n")
cata("<h3>Limma Analysis Output:</h3>\n")
cata("Links to PDF copies of plots are in 'Plots' section below <br />\n")

for (i in 1:nrow(imageData)) {
    if (grepl("density|box|mds|mdvol", imageData$Link[i])) {
        HtmlImage(imageData$Link[i], imageData$Label[i], width=1000)
    } else if (wantWeight) {
        HtmlImage(imageData$Link[i], imageData$Label[i], width=1000)
    } else {
        HtmlImage(imageData$Link[i], imageData$Label[i])
    }
}

cata("<h4>Differential Expression Counts:</h4>\n")

cata("<table border=\"1\" cellpadding=\"4\">\n")
cata("<tr>\n")
TableItem()
for (i in colnames(sigDiff)) {
    TableHeadItem(i)
}
cata("</tr>\n")
for (i in 1:nrow(sigDiff)) {
    cata("<tr>\n")
    TableHeadItem(unmake.names(row.names(sigDiff)[i]))
    for (j in 1:ncol(sigDiff)) {
        TableItem(as.character(sigDiff[i, j]))
    }
    cata("</tr>\n")
}
cata("</table>")

cata("<h4>Plots:</h4>\n")
#PDFs
for (i in 1:nrow(linkData)) {
    if (grepl("density|cpm|boxplot|mds|mdplots|voomplot|saplot", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

for (i in 1:nrow(linkData)) {
    if (grepl("mdplot_", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

for (i in 1:nrow(linkData)) {
    if (grepl("volplot", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

for (i in 1:nrow(linkData)) {
    if (grepl("heatmap", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

for (i in 1:nrow(linkData)) {
    if (grepl("stripcharts", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<h4>Tables:</h4>\n")
for (i in 1:nrow(linkData)) {
    if (grepl("counts$", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
    } else if (grepl(".tsv", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
    }
}

if (wantRda) {
    cata("<h4>R Data Object:</h4>\n")
    for (i in 1:nrow(linkData)) {
        if (grepl(".RData", linkData$Link[i])) {
            HtmlLink(linkData$Link[i], linkData$Label[i])
        }
    }
}

cata("<h4>Glimma Interactive Results:</h4>\n")
    for (i in 1:nrow(linkData)) {
        if (grepl("glimma", linkData$Link[i])) {
            HtmlLink(linkData$Link[i], linkData$Label[i])
        }
    }

cata("<p>Alt-click links to download file.</p>\n")
cata("<p>Click floppy disc icon associated history item to download ")
cata("all files.</p>\n")
cata("<p>.tsv files can be viewed in Excel or any spreadsheet program.</p>\n")

cata("<h4>Additional Information</h4>\n")
cata("<ul>\n")

if (filtCPM || filtSmpCount || filtTotCount) {
    if (filtCPM) {
    tempStr <- paste("Genes without more than", opt$cpmReq,
                                     "CPM in at least", opt$sampleReq, "samples are insignificant",
                                     "and filtered out.")
    } else if (filtSmpCount) {
        tempStr <- paste("Genes without more than", opt$cntReq,
                                     "counts in at least", opt$sampleReq, "samples are insignificant",
                                     "and filtered out.")
    } else if (filtTotCount) {
            tempStr <- paste("Genes without more than", opt$cntReq,
                                     "counts, after summing counts for all samples, are insignificant",
                                     "and filtered out.")
    }

    ListItem(tempStr)
    filterProp <- round(filteredCount/preFilterCount*100, digits=2)
    tempStr <- paste0(filteredCount, " of ", preFilterCount," (", filterProp,
                                     "%) genes were filtered out for low expression.")
    ListItem(tempStr)
}
ListItem(opt$normOpt, " was the method used to normalise library sizes.")
if (wantTrend) {
    ListItem("The limma-trend method was used.")
} else {
    ListItem("The limma-voom method was used.")
}
if (wantWeight) {
    ListItem("Weights were applied to samples.")
} else {
    ListItem("Weights were not applied to samples.")
}
if (wantTreat) {
    ListItem(paste("Testing significance relative to a fold-change threshold (TREAT) was performed using a threshold of log2 =", opt$lfcReq, "at FDR of", opt$pValReq, "."))
}
if (wantRobust) {
    if (wantTreat) {
        ListItem("TREAT was used with robust settings (robust=TRUE).")
    } else {
        ListItem("eBayes was used with robust settings (robust=TRUE).")
    }
}
if (opt$pAdjOpt!="none") {
    if (opt$pAdjOpt=="BH" || opt$pAdjOpt=="BY") {
        tempStr <- paste0("MD Plot highlighted genes are significant at FDR ",
                        "of ", opt$pValReq," and exhibit log2-fold-change of at ",
                        "least ", opt$lfcReq, ".")
        ListItem(tempStr)
    } else if (opt$pAdjOpt=="holm") {
        tempStr <- paste0("MD Plot highlighted genes are significant at adjusted ",
                        "p-value of ", opt$pValReq,"  by the Holm(1979) ",
                        "method, and exhibit log2-fold-change of at least ",
                        opt$lfcReq, ".")
        ListItem(tempStr)
    }
  } else {
        tempStr <- paste0("MD Plot highlighted genes are significant at p-value ",
                      "of ", opt$pValReq," and exhibit log2-fold-change of at ",
                      "least ", opt$lfcReq, ".")
        ListItem(tempStr)
}
cata("</ul>\n")

cata("<h4>Summary of experimental data:</h4>\n")

cata("<p>*CHECK THAT SAMPLES ARE ASSOCIATED WITH CORRECT GROUP(S)*</p>\n")

cata("<table border=\"1\" cellpadding=\"3\">\n")
cata("<tr>\n")
TableHeadItem("SampleID")
TableHeadItem(names(factors)[1]," (Primary Factor)")

if (ncol(factors) > 1) {
    for (i in names(factors)[2:length(names(factors))]) {
        TableHeadItem(i)
    }
    cata("</tr>\n")
}

for (i in 1:nrow(factors)) {
    cata("<tr>\n")
    TableHeadItem(row.names(factors)[i])
    for (j in 1:ncol(factors)) {
        TableItem(as.character(unmake.names(factors[i, j])))
  }
  cata("</tr>\n")
}
cata("</table>")

cit <- character()
link <- character()
link[1] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/limma/inst/doc/usersguide.pdf",
                  "\">", "limma User's Guide", "</a>.")

link[2] <- paste0("<a href=\"",
                  "http://www.bioconductor.org/packages/release/bioc/",
                  "vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf",
                  "\">", "edgeR User's Guide", "</a>")

cit[1] <- paste("Please cite the following paper for this tool:")

cit[2] <- paste("Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME,",
                "Asselin-Labat ML, Smyth GK, Ritchie ME (2015). Why weight? ",
                "Modelling sample and observational level variability improves power ",
                "in RNA-seq analyses. Nucleic Acids Research, 43(15), e97.")

cit[3] <- paste("Please cite the paper below for the limma software itself.",
                "Please also try to cite the appropriate methodology articles",
                "that describe the statistical methods implemented in limma,",
                "depending on which limma functions you are using. The",
                "methodology articles are listed in Section 2.1 of the",
                link[1],
                "Cite no. 3 only if sample weights were used.")
cit[4] <- paste("Smyth GK (2005). Limma: linear models for microarray data.",
                "In: 'Bioinformatics and Computational Biology Solutions using",
                "R and Bioconductor'. R. Gentleman, V. Carey, S. doit,.",
                "Irizarry, W. Huber (eds), Springer, New York, pages 397-420.")
cit[5] <- paste("Please cite the first paper for the software itself and the",
                "other papers for the various original statistical methods",
                "implemented in edgeR.  See Section 1.2 in the", link[2],
                "for more detail.")
cit[6] <- paste("Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a",
                "Bioconductor package for differential expression analysis",
                "of digital gene expression data. Bioinformatics 26, 139-140")
cit[7] <- paste("Robinson MD and Smyth GK (2007). Moderated statistical tests",
                "for assessing differences in tag abundance. Bioinformatics",
                "23, 2881-2887")
cit[8] <- paste("Robinson MD and Smyth GK (2008). Small-sample estimation of",
                "negative binomial dispersion, with applications to SAGE data.",
                "Biostatistics, 9, 321-332")
cit[9] <- paste("McCarthy DJ, Chen Y and Smyth GK (2012). Differential",
                "expression analysis of multifactor RNA-Seq experiments with",
                "respect to biological variation. Nucleic Acids Research 40,",
                "4288-4297")
cit[10] <- paste("Law CW, Chen Y, Shi W, and Smyth GK (2014). Voom:",
                "precision weights unlock linear model analysis tools for",
                "RNA-seq read counts. Genome Biology 15, R29.")
cit[11] <- paste("Ritchie ME, Diyagama D, Neilson J, van Laar R,",
                "Dobrovic A, Holloway A and Smyth GK (2006).",
                "Empirical array quality weights for microarray data.",
                "BMC Bioinformatics 7, Article 261.")
cata("<h3>Citations</h3>\n")
cata(cit[1], "\n")
cata("<br>\n")
cata(cit[2], "\n")

cata("<h4>limma</h4>\n")
cata(cit[3], "\n")
cata("<ol>\n")
ListItem(cit[4])
ListItem(cit[10])
ListItem(cit[11])
cata("</ol>\n")

cata("<h4>edgeR</h4>\n")
cata(cit[5], "\n")
cata("<ol>\n")
ListItem(cit[6])
ListItem(cit[7])
ListItem(cit[8])
ListItem(cit[9])
cata("</ol>\n")

cata("<p>Please report problems or suggestions to: su.s@wehi.edu.au</p>\n")

for (i in 1:nrow(linkData)) {
    if (grepl("session_info", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
    }
}

cata("<table border=\"0\">\n")
cata("<tr>\n")
TableItem("Task started at:"); TableItem(timeStart)
cata("</tr>\n")
cata("<tr>\n")
TableItem("Task ended at:"); TableItem(timeEnd)
cata("</tr>\n")
cata("<tr>\n")
TableItem("Task run time:"); TableItem(timeTaken)
cata("<tr>\n")
cata("</table>\n")

cata("</body>\n")
cata("</html>")
