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
#       lrtOpt", "t", 0, "logical"          -String specifying whether to perform LRT test instead 
#
# OUT: 
#       MDS Plot 
#       BCV Plot
#       QL Plot
#       MD Plot
#       Expression Table
#       HTML file linking to the ouputs
# Optional:
#       Normalised counts Table
#       RData file
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Oct 2017 (some code taken from the DESeq2 wrapper)

# Record starting time
timeStart <- as.character(Sys.time())

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Load all required libraries
library(methods, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)
library(scales, quietly=TRUE, warn.conflicts=FALSE)
library(getopt, quietly=TRUE, warn.conflicts=FALSE)

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
cata <- function(..., file=opt$htmlPath, sep="", fill=FALSE, labels=NULL, 
                                 append=TRUE) {
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
    "lrtOpt", "t", 0, "logical"),
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

if (is.null(opt$lrtOpt)) {
    wantLRT <- FALSE
} else {
    wantLRT <- TRUE
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
            # check samples names match
            if(!any(factorData[, 1] %in% colnames(counts)))
                stop("Sample IDs in factors file and count matrix don't match")
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

# Split up contrasts separated by comma into a vector then sanitise
contrastData <- unlist(strsplit(opt$contrastData, split=","))
contrastData <- sanitiseEquation(contrastData)
contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)

bcvOutPdf <- makeOut("bcvplot.pdf")
bcvOutPng <- makeOut("bcvplot.png")
qlOutPdf <- makeOut("qlplot.pdf")
qlOutPng <- makeOut("qlplot.png")
mdsOutPdf <- character()   # Initialise character vector
mdsOutPng <- character()
for (i in 1:ncol(factors)) {
    mdsOutPdf[i] <- makeOut(paste0("mdsplot_", names(factors)[i], ".pdf"))
    mdsOutPng[i] <- makeOut(paste0("mdsplot_", names(factors)[i], ".png"))
}
mdOutPdf <- character()
mdOutPng <- character()
topOut <- character()
for (i in 1:length(contrastData)) {
    mdOutPdf[i] <- makeOut(paste0("mdplot_", contrastData[i], ".pdf"))
    mdOutPng[i] <- makeOut(paste0("mdplot_", contrastData[i], ".png"))
    topOut[i] <- makeOut(paste0("edgeR_", contrastData[i], ".tsv"))
}   # Save output paths for each contrast as vectors
normOut <- makeOut("edgeR_normcounts.tsv")
rdaOut <- makeOut("edgeR_analysis.RData")
sessionOut <- makeOut("session_info.txt")

# Initialise data for html links and images, data frame with columns Label and 
# Link
linkData <- data.frame(Label=character(), Link=character(), stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(), stringsAsFactors=FALSE)

# Initialise vectors for storage of up/down/neutral regulated counts
upCount <- numeric()
downCount <- numeric()
flatCount <- numeric()

################################################################################
### Data Processing
################################################################################

# Extract counts and annotation data
data <- list()
data$counts <- counts
if (haveAnno) {
  # order annotation by genes in counts (assumes gene ids are in 1st column of geneanno)
  annoord <- geneanno[match(row.names(counts), geneanno[,1]), ]
  data$genes <- annoord
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
data$samples <- sampleanno
data$samples$lib.size <- colSums(data$counts)
data$samples$norm.factors <- 1
row.names(data$samples) <- colnames(data$counts)
data <- new("DGEList", data)

# Name rows of factors according to their sample
row.names(factors) <- names(data$counts)
factorList <- sapply(names(factors), pasteListName)

formula <- "~0" 
for (i in 1:length(factorList)) {
    formula <- paste(formula, factorList[i], sep="+")
}

formula <- formula(formula)
design <- model.matrix(formula)

for (i in 1:length(factorList)) {
    colnames(design) <- gsub(factorList[i], "", colnames(design), fixed=TRUE)
}

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method=opt$normOpt)

if (wantRobust) {
    data <- estimateDisp(data, design=design, robust=TRUE)
} else {
    data <- estimateDisp(data, design=design)
}

# Generate contrasts information
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

################################################################################
### Data Output
################################################################################

# Plot MDS
labels <- names(counts)

# MDS plot
png(mdsOutPng, width=600, height=600)
plotMDS(data, labels=labels, col=as.numeric(factors[, 1]), cex=0.8, main=paste("MDS Plot:", names(factors)[1]))
imgName <- paste0("MDS Plot_", names(factors)[1], ".png")
imgAddr <- paste0("mdsplot_", names(factors)[1], ".png")
imageData[1, ] <- c(imgName, imgAddr)
invisible(dev.off())

pdf(mdsOutPdf)
plotMDS(data, labels=labels, col=as.numeric(factors[, 1]), cex=0.8, main=paste("MDS Plot:", names(factors)[1]))
linkName <- paste0("MDS Plot_", names(factors)[1], ".pdf")
linkAddr <- paste0("mdsplot_", names(factors)[1], ".pdf")
linkData[1, ] <- c(linkName, linkAddr)
invisible(dev.off())

# If additional factors create additional MDS plots coloured by factor
if (ncol(factors) > 1) {
    for (i in 2:ncol(factors)) {
        png(mdsOutPng[i], width=600, height=600)
        plotMDS(data, labels=labels, col=as.numeric(factors[, i]), cex=0.8, main=paste("MDS Plot:", names(factors)[i]))
        imgName <- paste0("MDS Plot_", names(factors)[i], ".png")
        imgAddr <- paste0("mdsplot_", names(factors)[i], ".png")
        imageData <- rbind(imageData, c(imgName, imgAddr))
        invisible(dev.off())

        pdf(mdsOutPdf[i])
        plotMDS(data, labels=labels, col=as.numeric(factors[, i]), cex=0.8, main=paste("MDS Plot:", names(factors)[i]))
        linkName <- paste0("MDS Plot_", names(factors)[i], ".pdf")
        linkAddr <- paste0("mdsplot_", names(factors)[i], ".pdf")
        linkData <- rbind(linkData, c(linkName, linkAddr))
        invisible(dev.off())
    }
}

# BCV Plot
png(bcvOutPng, width=600, height=600)
plotBCV(data, main="BCV Plot")
imgName <- "BCV Plot"
imgAddr <- "bcvplot.png"
imageData <- rbind(imageData, c(imgName, imgAddr))
invisible(dev.off())

pdf(bcvOutPdf)
plotBCV(data, main="BCV Plot")
linkName <- paste0("BCV Plot.pdf")
linkAddr <- paste0("bcvplot.pdf")
linkData <- rbind(linkData, c(linkName, linkAddr))
invisible(dev.off())

# Generate fit
if (wantLRT) {
    
    fit <- glmFit(data, design)
    
} else {
    
    if (wantRobust) {
        fit <- glmQLFit(data, design, robust=TRUE)
    } else {
        fit <- glmQLFit(data, design)
    }

    # Plot QL dispersions
    png(qlOutPng, width=600, height=600)
    plotQLDisp(fit, main="QL Plot")
    imgName <- "QL Plot"
    imgAddr <- "qlplot.png"
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())

    pdf(qlOutPdf)
    plotQLDisp(fit, main="QL Plot")
    linkName <- "QL Plot.pdf"
    linkAddr <- "qlplot.pdf"
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())
}

 # Save normalised counts (log2cpm)
if (wantNorm) { 
        normalisedCounts <- cpm(data, normalized.lib.sizes=TRUE, log=TRUE) 
        normalisedCounts <- data.frame(data$genes, normalisedCounts)
        write.table (normalisedCounts, file=normOut, row.names=FALSE, sep="\t", quote=FALSE)
        linkData <- rbind(linkData, c("edgeR_normcounts.tsv", "edgeR_normcounts.tsv"))
}


for (i in 1:length(contrastData)) {
    if (wantLRT) {
        res <- glmLRT(fit, contrast=contrasts[, i])
    } else {
        res <- glmQLFTest(fit, contrast=contrasts[, i])
    }

    status = decideTestsDGE(res, adjust.method=opt$pAdjOpt, p.value=opt$pValReq,
                                             lfc=opt$lfcReq)
    sumStatus <- summary(status)

    # Collect counts for differential expression
    upCount[i] <- sumStatus["Up", ]
    downCount[i] <- sumStatus["Down", ]
    flatCount[i] <- sumStatus["NotSig", ]
                                             
    # Write top expressions table
    top <- topTags(res, adjust.method=opt$pAdjOpt, n=Inf, sort.by="PValue")
    write.table(top, file=topOut[i], row.names=FALSE, sep="\t", quote=FALSE)
    
    linkName <- paste0("edgeR_", contrastData[i], ".tsv")
    linkAddr <- paste0("edgeR_", contrastData[i], ".tsv")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    
    # Plot MD (log ratios vs mean difference) using limma package
    pdf(mdOutPdf[i])
    limma::plotMD(res, status=status,
                                main=paste("MD Plot:", unmake.names(contrastData[i])), 
                                hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1),
                                xlab="Average Expression", ylab="logFC")
    
    abline(h=0, col="grey", lty=2)
    
    linkName <- paste0("MD Plot_", contrastData[i], ".pdf")
    linkAddr <- paste0("mdplot_", contrastData[i], ".pdf")
    linkData <- rbind(linkData, c(linkName, linkAddr))
    invisible(dev.off())
    
    png(mdOutPng[i], height=600, width=600)
    limma::plotMD(res, status=status,
                                main=paste("MD Plot:", unmake.names(contrastData[i])), 
                                hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1),
                                xlab="Average Expression", ylab="logFC")
    
    abline(h=0, col="grey", lty=2)
    
    imgName <- paste0("MD Plot_", contrastData[i], ".png")
    imgAddr <- paste0("mdplot_", contrastData[i], ".png")
    imageData <- rbind(imageData, c(imgName, imgAddr))
    invisible(dev.off())
}
sigDiff <- data.frame(Up=upCount, Flat=flatCount, Down=downCount)
row.names(sigDiff) <- contrastData

# Save relevant items as rda object
if (wantRda) {
    if (wantNorm) {
        save(counts, data, status, normalisedCounts, labels, factors, fit, res, top, contrasts, design,
                 file=rdaOut, ascii=TRUE)
    } else {
        save(counts, data, status, labels, factors, fit, res, top, contrasts, design,
                 file=rdaOut, ascii=TRUE)
    }
    linkData <- rbind(linkData, c("edgeR_analysis.RData", "edgeR_analysis.RData"))
}

# Record session info
writeLines(capture.output(sessionInfo()), sessionOut)
linkData <- rbind(linkData, c("Session Info", "session_info.txt"))

# Record ending time and calculate total run time
timeEnd <- as.character(Sys.time())
timeTaken <- capture.output(round(difftime(timeEnd, timeStart), digits=3))
timeTaken <- gsub("Time difference of ", "", timeTaken, fixed=TRUE)

################################################################################
### HTML Generation
################################################################################

# Clear file
cat("", file=opt$htmlPath)

cata("<html>\n")

cata("<body>\n")
cata("<h3>edgeR Analysis Output:</h3>\n")
cata("Links to PDF copies of plots are in 'Plots' section below.<br />\n")

HtmlImage(imageData$Link[1], imageData$Label[1])

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
    cata("<h4>R Data Objects:</h4>\n")
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
if (wantLRT) {
    ListItem("The edgeR likelihood ratio test was used.")
} else {
    if (wantRobust) {
        ListItem("The edgeR quasi-likelihood test was used with robust settings (robust=TRUE with estimateDisp and glmQLFit).")
    } else {
            ListItem("The edgeR quasi-likelihood test was used.")
    }
}
if (opt$pAdjOpt!="none") {
    if (opt$pAdjOpt=="BH" || opt$pAdjOpt=="BY") {
        tempStr <- paste0("MD-Plot highlighted genes are significant at FDR ",
                                            "of ", opt$pValReq," and exhibit log2-fold-change of at ", 
                                            "least ", opt$lfcReq, ".")
        ListItem(tempStr)
    } else if (opt$pAdjOpt=="holm") {
        tempStr <- paste0("MD-Plot highlighted genes are significant at adjusted ",
                                            "p-value of ", opt$pValReq,"  by the Holm(1979) ",
                                            "method, and exhibit log2-fold-change of at least ", 
                                            opt$lfcReq, ".")
        ListItem(tempStr)
    }
} else {
    tempStr <- paste0("MD-Plot highlighted genes are significant at p-value ",
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
TableHeadItem(names(factors)[1], " (Primary Factor)")

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
