# This tool takes in a matrix of feature counts as well as gene annotations and
# outputs a table of top expressions as well as various plots for differential
# expression analysis
#
# ARGS: 1.countPath       -Path to RData input containing counts
#       2.annoPath        -Path to input containing gene annotations
#       3.htmlPath        -Path to html file linking to other outputs
#       4.outPath         -Path to folder to write all output to
#       5.rdaOpt          -String specifying if RData should be saved
#       6.normOpt         -String specifying type of normalisation used
#       7.weightOpt       -String specifying usage of weights
#       8.contrastData    -String containing contrasts of interest
#       9.cpmReq          -Float specifying cpm requirement
#       10.sampleReq      -Integer specifying cpm requirement
#       11.pAdjOpt        -String specifying the p-value adjustment method
#       12.pValReq        -Float specifying the p-value requirement
#       13.lfcReq         -Float specifying the log-fold-change requirement
#       14.normCounts     -String specifying if normalised counts should be output
#       15.factPath       -Path to factor information file
#       16.factorData     -Strings containing factor names and values if manually input 
#
# OUT:  Voom Plot
#       BCV Plot
#       MA Plot
#       Expression Table
#       HTML file linking to the ouputs
#
# Author: Shian Su - registertonysu@gmail.com - Jan 2014
# Modified by: Maria Doyle - Jun 2017

# Record starting time
timeStart <- as.character(Sys.time())

# Load all required libraries
library(methods, quietly=TRUE, warn.conflicts=FALSE)
library(statmod, quietly=TRUE, warn.conflicts=FALSE)
library(splines, quietly=TRUE, warn.conflicts=FALSE)
library(edgeR, quietly=TRUE, warn.conflicts=FALSE)
library(limma, quietly=TRUE, warn.conflicts=FALSE)
library(scales, quietly=TRUE, warn.conflicts=FALSE)

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
  return(paste0(outPath, "/", filename))
}

# Generating design information
pasteListName <- function(string) {
  return(paste0("factors$", string))
}

# Create cata function: default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)
cata <- function(..., file = htmlPath, sep = "", fill = FALSE, labels = NULL, 
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

# Collects arguments from command line
argv <- commandArgs(TRUE)

# Grab arguments
countPath <- as.character(argv[1])
annoPath <- as.character(argv[2])
htmlPath <- as.character(argv[3])
outPath <- as.character(argv[4])
rdaOpt <- as.character(argv[5])
normOpt <- as.character(argv[6])
weightOpt <- as.character(argv[7])
contrastData <- as.character(argv[8])
cpmReq <- as.numeric(argv[9])
sampleReq <- as.numeric(argv[10])
pAdjOpt <- as.character(argv[11])
pValReq <- as.numeric(argv[12])
lfcReq <- as.numeric(argv[13])
normCounts <- as.character(argv[14])
factPath <- as.character(argv[15])
# Process factors
if (as.character(argv[16])=="None") {
    factorData <- read.table(factPath, header=TRUE, sep="\t")
    factors <- factorData[,-1, drop=FALSE]
}  else { 
    factorData <- list()
    for (i in 16:length(argv)) {
        newFact <- unlist(strsplit(as.character(argv[i]), split="::"))
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

# Process other arguments
if (weightOpt=="yes") {
  wantWeight <- TRUE
} else {
  wantWeight <- FALSE
}

if (rdaOpt=="yes") {
  wantRda <- TRUE
} else {
  wantRda <- FALSE
}

if (annoPath=="None") {
  haveAnno <- FALSE
} else {
  haveAnno <- TRUE
}

if (normCounts=="yes") {
  wantNorm <- TRUE
} else {
  wantNorm <- FALSE
}


#Create output directory
dir.create(outPath, showWarnings=FALSE)

# Split up contrasts seperated by comma into a vector then sanitise
contrastData <- unlist(strsplit(contrastData, split=","))
contrastData <- sanitiseEquation(contrastData)
contrastData <- gsub(" ", ".", contrastData, fixed=TRUE)

bcvOutPdf <- makeOut("bcvplot.pdf")
bcvOutPng <- makeOut("bcvplot.png")
mdsOutPdf <- makeOut("mdsplot.pdf")
mdsOutPng <- makeOut("mdsplot.png")
voomOutPdf <- makeOut("voomplot.pdf")
voomOutPng <- makeOut("voomplot.png") 
maOutPdf <- character()   # Initialise character vector
maOutPng <- character()
topOut <- character()
for (i in 1:length(contrastData)) {
  maOutPdf[i] <- makeOut(paste0("maplot_", contrastData[i], ".pdf"))
  maOutPng[i] <- makeOut(paste0("maplot_", contrastData[i], ".png"))
  topOut[i] <- makeOut(paste0("limma-voom_", contrastData[i], ".tsv"))
}                         # Save output paths for each contrast as vectors
normOut <- makeOut("limma-voom_normcounts.tsv")
rdaOut <- makeOut("RData.rda")
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
                        
# Read in counts and geneanno data
counts <- read.table(countPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
row.names(counts) <- counts[, 1]
counts <- counts[ , -1]
countsRows <- nrow(counts)
if (haveAnno) {
  geneanno <- read.table(annoPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
}

################################################################################
### Data Processing
################################################################################

# Extract counts and annotation data
data <- list()
data$counts <- counts
if (haveAnno) {
  data$genes <- geneanno
} else {
  data$genes <- data.frame(GeneID=row.names(counts))
}

# Filter out genes that do not have a required cpm in a required number of
# samples
preFilterCount <- nrow(data$counts)
sel <- rowSums(cpm(data$counts) > cpmReq) >= sampleReq
data$counts <- data$counts[sel, ]
data$genes <- data$genes[sel, ,drop = FALSE]
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

# Calculating normalising factor, estimating dispersion
data <- calcNormFactors(data, method=normOpt)
#data <- estimateDisp(data, design=design, robust=TRUE)

# Generate contrasts information
contrasts <- makeContrasts(contrasts=contrastData, levels=design)

# Name rows of factors according to their sample
row.names(factors) <- names(data$counts)

################################################################################
### Data Output
################################################################################

# BCV Plot
#png(bcvOutPng, width=600, height=600)
#plotBCV(data, main="BCV Plot")
#imageData[1, ] <- c("BCV Plot", "bcvplot.png")
#invisible(dev.off())

#pdf(bcvOutPdf)
#plotBCV(data, main="BCV Plot")
#invisible(dev.off())

if (wantWeight) {
  # Creating voom data object and plot
  png(voomOutPng, width=1000, height=600)
  vData <- voomWithQualityWeights(data, design=design, plot=TRUE)
  imageData[1, ] <- c("Voom Plot", "voomplot.png")
  invisible(dev.off())
  
  pdf(voomOutPdf, width=14)
  vData <- voomWithQualityWeights(data, design=design, plot=TRUE)
  linkData[1, ] <- c("Voom Plot (.pdf)", "voomplot.pdf")
  invisible(dev.off())
  
  # Generating fit data and top table with weights
  wts <- vData$weights
  voomFit <- lmFit(vData, design, weights=wts)
  
} else {
  # Creating voom data object and plot
  png(voomOutPng, width=600, height=600)
  vData <- voom(data, design=design, plot=TRUE)
  imageData[1, ] <- c("Voom Plot", "voomplot.png")
  invisible(dev.off())
  
  pdf(voomOutPdf)
  vData <- voom(data, design=design, plot=TRUE)
  linkData[1, ] <- c("Voom Plot (.pdf)", "voomplot.pdf")
  invisible(dev.off())
  
  # Generate voom fit
  voomFit <- lmFit(vData, design)
  
}

 # Save normalised counts (log2cpm)
if (wantNorm) {   
    norm_counts <- data.frame(vData$genes, vData$E)
    write.table (norm_counts, file=normOut, row.names=FALSE, sep="\t")
    linkData <- rbind(linkData, c("limma-voom_normcounts.tsv", "limma-voom_normcounts.tsv"))
}

# Fit linear model and estimate dispersion with eBayes
voomFit <- contrasts.fit(voomFit, contrasts)
voomFit <- eBayes(voomFit)

# Plot MDS
labels <- names(counts)
png(mdsOutPng, width=600, height=600)
# Currently only using a single factor
plotMDS(vData, labels=labels, col=as.numeric(factors[, 1]), cex=0.8)
imgName <- "Voom Plot"
imgAddr <- "mdsplot.png"
imageData <- rbind(imageData, c(imgName, imgAddr))
invisible(dev.off())

pdf(mdsOutPdf)
plotMDS(vData, labels=labels, cex=0.5)
linkName <- paste0("MDS Plot (.pdf)")
linkAddr <- paste0("mdsplot.pdf")
linkData <- rbind(linkData, c(linkName, linkAddr))
invisible(dev.off())


for (i in 1:length(contrastData)) {

  status = decideTests(voomFit[, i], adjust.method=pAdjOpt, p.value=pValReq,
                       lfc=lfcReq)
                       
  sumStatus <- summary(status)
  
  # Collect counts for differential expression
  upCount[i] <- sumStatus["1",]
  downCount[i] <- sumStatus["-1",]
  flatCount[i] <- sumStatus["0",]
                       
  # Write top expressions table
  top <- topTable(voomFit, coef=i, number=Inf, sort.by="P")
  write.table(top, file=topOut[i], row.names=FALSE, sep="\t")
  
  linkName <- paste0("limma-voom_", contrastData[i], ".tsv")
  linkAddr <- paste0("limma-voom_", contrastData[i], ".tsv")
  linkData <- rbind(linkData, c(linkName, linkAddr))
  
  # Plot MA (log ratios vs mean average) using limma package on weighted 
  pdf(maOutPdf[i])
  limma::plotMA(voomFit, status=status, coef=i,
                main=paste("MA Plot:", unmake.names(contrastData[i])), 
                col=alpha(c("firebrick", "blue"), 0.4), values=c("1", "-1"),
                xlab="Average Expression", ylab="logFC")
  
  abline(h=0, col="grey", lty=2)
  
  linkName <- paste0("MA Plot_", contrastData[i], " (.pdf)")
  linkAddr <- paste0("maplot_", contrastData[i], ".pdf")
  linkData <- rbind(linkData, c(linkName, linkAddr))
  invisible(dev.off())
  
  png(maOutPng[i], height=600, width=600)
  limma::plotMA(voomFit, status=status, coef=i,
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
  if (wantWeight) {
    save(data, status, vData, labels, factors, wts, voomFit, top, contrasts, 
         design,
         file=rdaOut, ascii=TRUE)
  } else {
    save(data, status, vData, labels, factors, voomFit, top, contrasts, design,
         file=rdaOut, ascii=TRUE)
  }
  linkData <- rbind(linkData, c("RData (.rda)", "RData.rda"))
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
cat("", file=htmlPath)

cata("<html>\n")

cata("<body>\n")
cata("<h3>Limma-voom Analysis Output:</h3>\n")
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
    if (grepl(".rda", linkData$Link[i])) {
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
if (cpmReq!=0 && sampleReq!=0) {
  tempStr <- paste("Genes without more than", cpmReq,
                   "CPM in at least", sampleReq, "samples are insignificant",
                   "and filtered out.")
  ListItem(tempStr)
  filterProp <- round(filteredCount/preFilterCount*100, digits=2)
  tempStr <- paste0(filteredCount, " of ", preFilterCount," (", filterProp,
                   "%) genes were filtered out for low expression.")
  ListItem(tempStr)
}
ListItem(normOpt, " was the method used to normalise library sizes.")
if (wantWeight) {
  ListItem("Weights were applied to samples.")
} else {
  ListItem("Weights were not applied to samples.")
}
if (pAdjOpt!="none") {
  if (pAdjOpt=="BH" || pAdjOpt=="BY") {
    tempStr <- paste0("MA-Plot highlighted genes are significant at FDR ",
                      "of ", pValReq," and exhibit log2-fold-change of at ", 
                      "least ", lfcReq, ".")
    ListItem(tempStr)
  } else if (pAdjOpt=="holm") {
    tempStr <- paste0("MA-Plot highlighted genes are significant at adjusted ",
                      "p-value of ", pValReq,"  by the Holm(1979) ",
                      "method, and exhibit log2-fold-change of at least ", 
                      lfcReq, ".")
    ListItem(tempStr)
  }
} else {
  tempStr <- paste0("MA-Plot highlighted genes are significant at p-value ",
                    "of ", pValReq," and exhibit log2-fold-change of at ", 
                    "least ", lfcReq, ".")
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
