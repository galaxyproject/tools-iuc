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
#
# OUT:
#       MDS Plot
#       Voom/SA plot
#       MD Plot
#       Expression Table
#       HTML file linking to the ouputs
# Optional:
#       Normalised counts Table
#       RData file
#
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Jun 2017, Jan 2018

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

if (packageVersion("limma") < "3.20.1") {
    stop("Please update 'limma' to version >= 3.20.1 to run this tool")
}

################################################################################
### Function Delcaration
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
HtmlImage <- function(source, label=source, height=600, width=600) {
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
    "normCounts", "x", 0, "logical",
    "rdaOpt", "r", 0, "logical",
    "lfcReq", "l", 1, "double",
    "pValReq", "p", 1, "double",
    "pAdjOpt", "d", 1, "character",
    "normOpt", "n", 1, "character",
    "robOpt", "b", 0, "logical",
    "trend", "t", 1, "double",
    "weightOpt", "w", 0, "logical"),
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
    counts <- read.table(opt$matrixPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    row.names(counts) <- counts[, 1]
    counts <- counts[ , -1]
    countsRows <- nrow(counts)

    # Process factors
    if (is.null(opt$factInput)) {
            factorData <- read.table(opt$factFile, header=TRUE, sep="\t")
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
    geneanno <- read.table(opt$annoPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}

#Create output directory
dir.create(opt$outPath, showWarnings=FALSE)

# Split up contrasts seperated by comma into a vector then sanitise
contrastData <- unlist(strsplit(opt$contrastData, split=","))
contrastData <- sanitiseEquation(contrastData)
contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)


mdsOutPdf <- makeOut("mdsplot_nonorm.pdf")
mdsOutPng <- makeOut("mdsplot_nonorm.png")
nmdsOutPdf <- makeOut("mdsplot.pdf")
nmdsOutPng <- makeOut("mdsplot.png")
maOutPdf <- character()   # Initialise character vector
maOutPng <- character()
topOut <- character()
for (i in 1:length(contrastData)) {
    maOutPdf[i] <- makeOut(paste0("maplot_", contrastData[i], ".pdf"))
    maOutPng[i] <- makeOut(paste0("maplot_", contrastData[i], ".png"))
    topOut[i] <- makeOut(paste0(deMethod, "_", contrastData[i], ".tsv"))
}
normOut <- makeOut(paste0(deMethod, "_normcounts.tsv"))
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
  data$genes <- geneanno
} else {
  data$genes <- data.frame(GeneID=row.names(counts))
}

# If filter crieteria set, filter out genes that do not have a required cpm/counts in a required number of
# samples. Default is no filtering
preFilterCount <- nrow(data$counts)

if (filtCPM || filtSmpCount || filtTotCount) {

    if (filtTotCount) {
        keep <- rowSums(data$counts) >= opt$cntReq
    } else if (filtSmpCount) {
        keep <- rowSums(data$counts >= opt$cntReq) >= opt$sampleReq
    } else if (filtCPM) {
        keep <- rowSums(cpm(data$counts) >= opt$cpmReq) >= opt$sampleReq
    }

    data$counts <- data$counts[keep, ]
    data$genes <- data$genes[keep, , drop=FALSE]
}

postFilterCount <- nrow(data$counts)
filteredCount <- preFilterCount-postFilterCount

# Creating naming data
samplenames <- colnames(data$counts)
sampleanno <- data.frame("sampleID"=samplenames, factors)


# Generating the DGEList object "data"
print("Generating DGEList object")
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

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
data <- calcNormFactors(data, method=opt$normOpt)

# Generate contrasts information
print("Generating Contrasts")
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

################################################################################
### Data Output
################################################################################
# Plot MDS
print("Generating MDS plot")
labels <- names(counts)
png(mdsOutPng, width=600, height=600)
# Currently only using a single factor
plotMDS(data, labels=labels, col=as.numeric(factors[, 1]), cex=0.8, main="MDS Plot (unnormalised)")
imageData[1, ] <- c("MDS Plot (unnormalised)", "mdsplot_nonorm.png")
invisible(dev.off())

pdf(mdsOutPdf)
plotMDS(data, labels=labels, cex=0.5)
linkData[1, ] <- c("MDS Plot (unnormalised).pdf", "mdsplot_nonorm.pdf")
invisible(dev.off())

if (wantTrend) {
    # limma-trend approach
    logCPM <- cpm(data, log=TRUE, prior.count=opt$trend)
    fit <- lmFit(logCPM, design)
    fit <- contrasts.fit(fit, contrasts)
    if (wantRobust) {
        fit <- eBayes(fit, trend=TRUE, robust=TRUE)
    } else {
        fit <- eBayes(fit, trend=TRUE, robust=FALSE)
    }
    # plot fit with plotSA
    saOutPng <- makeOut("saplot.png")
    saOutPdf <- makeOut("saplot.pdf")

    png(saOutPng, width=600, height=600)
    plotSA(fit, main="SA Plot")
    imgName <- "SA Plot.png"
    imgAddr <- "saplot.png"
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())

    pdf(saOutPdf, width=14)
    plotSA(fit, main="SA Plot")
    linkName <- paste0("SA Plot.pdf")
    linkAddr <- paste0("saplot.pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())

    plotData <- logCPM

    # Save normalised counts (log2cpm)
    if (wantNorm) {
        write.table(logCPM, file=normOut, row.names=TRUE, sep="\t")
        linkData <- rbind(linkData, c((paste0(deMethod, "_", "normcounts.tsv")), (paste0(deMethod, "_", "normcounts.tsv"))))
    }
} else {
    # limma-voom approach
    voomOutPdf <- makeOut("voomplot.pdf")
    voomOutPng <- makeOut("voomplot.png")

    if (wantWeight) {
        # Creating voom data object and plot
        png(voomOutPng, width=1000, height=600)
        vData <- voomWithQualityWeights(data, design=design, plot=TRUE)
        imgName <- "Voom Plot.png"
        imgAddr <- "voomplot.png"
        imageData <- rbind(imageData, c(imgName, imgAddr))
        invisible(dev.off())

        pdf(voomOutPdf, width=14)
        vData <- voomWithQualityWeights(data, design=design, plot=TRUE)
        linkName <- paste0("Voom Plot.pdf")
        linkAddr <- paste0("voomplot.pdf")
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())

        # Generating fit data and top table with weights
        wts <- vData$weights
        voomFit <- lmFit(vData, design, weights=wts)

    } else {
        # Creating voom data object and plot
        png(voomOutPng, width=600, height=600)
        vData <- voom(data, design=design, plot=TRUE)
        imgName <- "Voom Plot"
        imgAddr <- "voomplot.png"
        imageData <- rbind(imageData, c(imgName, imgAddr))
        invisible(dev.off())

        pdf(voomOutPdf)
        vData <- voom(data, design=design, plot=TRUE)
        linkName <- paste0("Voom Plot.pdf")
        linkAddr <- paste0("voomplot.pdf")
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())

        # Generate voom fit
        voomFit <- lmFit(vData, design)
    }

     # Save normalised counts (log2cpm)
    if (wantNorm) {
        norm_counts <- data.frame(vData$genes, vData$E)
        write.table(norm_counts, file=normOut, row.names=FALSE, sep="\t")
        linkData <- rbind(linkData, c((paste0(deMethod, "_", "normcounts.tsv")), (paste0(deMethod, "_", "normcounts.tsv"))))
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

print("Generating normalised MDS plot")
png(nmdsOutPng, width=600, height=600)
# Currently only using a single factor
plotMDS(plotData, labels=labels, col=as.numeric(factors[, 1]), cex=0.8, main="MDS Plot (normalised)")
imgName <- "MDS Plot (normalised)"
imgAddr <- "mdsplot.png"
imageData <- rbind(imageData, c(imgName, imgAddr))
invisible(dev.off())

pdf(nmdsOutPdf)
plotMDS(plotData, labels=labels, cex=0.5)
linkName <- paste0("MDS Plot (normalised).pdf")
linkAddr <- paste0("mdsplot.pdf")
linkData <- rbind(linkData, c(linkName, linkAddr))
invisible(dev.off())


print("Generating DE results")
status = decideTests(fit, adjust.method=opt$pAdjOpt, p.value=opt$pValReq,
                       lfc=opt$lfcReq)
sumStatus <- summary(status)

for (i in 1:length(contrastData)) {
    # Collect counts for differential expression
    upCount[i] <- sumStatus["Up", i]
    downCount[i] <- sumStatus["Down", i]
    flatCount[i] <- sumStatus["NotSig", i]

    # Write top expressions table
    top <- topTable(fit, coef=i, number=Inf, sort.by="P")
    if (wantTrend) {
        write.table(top, file=topOut[i], row.names=TRUE, sep="\t")
    } else {
        write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
    }

    linkName <- paste0(deMethod, "_", contrastData[i], ".tsv")
    linkAddr <- paste0(deMethod, "_", contrastData[i], ".tsv")
    linkData <- rbind(linkData, c(linkName, linkAddr))

    # Plot MA (log ratios vs mean average) using limma package on weighted
    pdf(maOutPdf[i])
    limma::plotMD(fit, status=status, coef=i,
                  main=paste("MA Plot:", unmake.names(contrastData[i])),
                  col=alpha(c("firebrick", "blue"), 0.4), values=c("1", "-1"),
                  xlab="Average Expression", ylab="logFC")

    abline(h=0, col="grey", lty=2)

    linkName <- paste0("MA Plot_", contrastData[i], " (.pdf)")
    linkAddr <- paste0("maplot_", contrastData[i], ".pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())

    png(maOutPng[i], height=600, width=600)
    limma::plotMD(fit, status=status, coef=i,
                  main=paste("MA Plot:", unmake.names(contrastData[i])),
                  col=alpha(c("firebrick", "blue"), 0.4), values=c("1", "-1"),
                  xlab="Average Expression", ylab="logFC")

    abline(h=0, col="grey", lty=2)

    imgName <- paste0("MA Plot_", contrastData[i])
    imgAddr <- paste0("maplot_", contrastData[i], ".png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())
}
sigDiff <- data.frame(Up=upCount, Flat=flatCount, Down=downCount)
row.names(sigDiff) <- contrastData

# Save relevant items as rda object
if (wantRda) {
    print("Saving RData")
    if (wantWeight) {
      save(data, status, plotData, labels, factors, wts, fit, top, contrasts,
           design,
           file=rdaOut, ascii=TRUE)
    } else {
      save(data, status, plotData, labels, factors, fit, top, contrasts, design,
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
cata("PDF copies of JPEGS available in 'Plots' section.<br />\n")
if (wantWeight) {
    HtmlImage(imageData$Link[1], imageData$Label[1], width=1000)
} else {
    HtmlImage(imageData$Link[1], imageData$Label[1])
}

for (i in 2:nrow(imageData)) {
    HtmlImage(imageData$Link[i], imageData$Label[i])
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
for (i in 1:nrow(linkData)) {
    if (grepl(".pdf", linkData$Link[i])) {
        HtmlLink(linkData$Link[i], linkData$Label[i])
  }
}

cata("<h4>Tables:</h4>\n")
for (i in 1:nrow(linkData)) {
    if (grepl(".tsv", linkData$Link[i])) {
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

cata("<p>Alt-click links to download file.</p>\n")
cata("<p>Click floppy disc icon associated history item to download ")
cata("all files.</p>\n")
cata("<p>.tsv files can be viewed in Excel or any spreadsheet program.</p>\n")

cata("<h4>Additional Information</h4>\n")
cata("<ul>\n")

if (filtCPM || filtSmpCount || filtTotCount) {
    if (filtCPM) {
    tempStr <- paste("Genes without more than", opt$cmpReq,
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
if (wantRobust) {
    ListItem("eBayes was used with robust settings (robust=TRUE).")
}
if (opt$pAdjOpt!="none") {
    if (opt$pAdjOpt=="BH" || opt$pAdjOpt=="BY") {
        tempStr <- paste0("MA-Plot highlighted genes are significant at FDR ",
                        "of ", opt$pValReq," and exhibit log2-fold-change of at ",
                        "least ", opt$lfcReq, ".")
        ListItem(tempStr)
    } else if (opt$pAdjOpt=="holm") {
        tempStr <- paste0("MA-Plot highlighted genes are significant at adjusted ",
                        "p-value of ", opt$pValReq,"  by the Holm(1979) ",
                        "method, and exhibit log2-fold-change of at least ",
                        opt$lfcReq, ".")
        ListItem(tempStr)
    }
  } else {
        tempStr <- paste0("MA-Plot highlighted genes are significant at p-value ",
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
