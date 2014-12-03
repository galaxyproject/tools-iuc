## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library('getopt')
library('rjson')
library('DESeq2')
library("RColorBrewer")
library("gplots")

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help' , 'h', 0, "logical",
    'outfile' , 'o', 1, "character",
    'outfilefiltered' , 'f', 1, "character",
    'plots' , 'p', 2, "character",
    'factors', 'm', 2, "character",
    'filtermode', 'l', 2, "character",
    'threshold', 'c', 2, "double",
    'organism', 'g', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

if( is.null(opt$filtermode))
    opt$filtermode = "absolute"

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)


generatePlots <-function(dds, title_prefix){
    plotDispEsts(dds, main= paste(title_prefix, "Dispersion estimate plot"))
    plotMA(dds, main= paste(title_prefix, "MA-plot"))

    select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    rld <- rlogTransformation(dds, blind=TRUE)
    vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    heatmap.2(mat, trace="none", col = rev(hmcol), main = paste(title_prefix, "Sample-to-sample distances"), margin=c(13,13))
}

parser <- newJSONParser()
parser$addData( opt$factors )
factorsList <- parser$getObject()

sampleTable<-data.frame()
factorNames<-c()
primaryfactor = TRUE
for(factor in factorsList){
    factorName<-factor[[1]]
    factorNames<-append(factorNames, factorName)
    factorValuesMapList<-factor[[2]]
    c = length(factorValuesMapList)
    for (factorValuesMap in factorValuesMapList){
        for(files in factorValuesMap){
            fvc = 0
            for(file in files){
                fvc = fvc+1
                if(primaryfactor) {
                    sampleTable[basename(file),"sampleName"]<-paste(fvc,names(factorValuesMap),sep="_")
                }
                sampleTable[basename(file),"fileName"]<-file
                sampleTable[basename(file),factorName]<-paste(c,names(factorValuesMap),sep="_")
            }
        }
        c = c-1
    }
    primaryfactor = FALSE
}

factorNames<-rev(factorNames)
designFormula <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = "",
                                      design =  designFormula)

ddsHTSeq <- DESeq(ddsHTSeq)
#ddsHTSeq <- nbinomWaldTest(ddsHTSeq, cooksCutoff=FALSE)
deseqRes <- results(ddsHTSeq)
mcols(deseqRes)$description

resSorted <- deseqRes[order(deseqRes$padj),]
head(resSorted)
write.table(as.data.frame(resSorted), file = opt$outfile, sep="\t", quote = FALSE, append=TRUE, col.names = FALSE)


if(opt$filtermode == "absolute"){
    use = (deseqRes$baseMean > opt$threshold & !is.na(deseqRes$pvalue))
} else if(opt$filtermode == "quantile"){
    opt$threshold = opt$threshold/100
    print(paste("Threshold:", opt$threshold))
    print(quantile(unique(deseqRes$baseMean), probs = opt$threshold))
    use = (deseqRes$baseMean > quantile(unique(deseqRes$baseMean), probs = opt$threshold) & (!is.na(deseqRes$pvalue)))
} else{
    print(paste("Unknown filtermode:", opt$filtermode))
    q("no", 1, FALSE)
}

ddsHTSeqFilt <- ddsHTSeq[use, ]
ddsHTSeqFilt <- DESeq(ddsHTSeqFilt)
#ddsHTSeqFilt <- nbinomWaldTest(ddsHTSeqFilt, cooksCutoff=FALSE)
deseqResFilt <- results(ddsHTSeqFilt)
mcols(deseqResFilt)$description

resFiltSorted <- deseqResFilt[order(deseqResFilt$padj),]
head(resFiltSorted)
write.table(as.data.frame(resFiltSorted), file = opt$outfilefiltered, sep="\t", quote = FALSE, append=TRUE, col.names = FALSE)

if(!is.null(opt$plots)){
    pdf(opt$plots)
    generatePlots(ddsHTSeq, "Complete: ")
    generatePlots(ddsHTSeqFilt, "Filtered: ")
    ## p-value distribution
    h1 <- hist(deseqRes$pvalue[!use], breaks=50, plot=FALSE)
    h2 <- hist(deseqRes$pvalue[use], breaks=50, plot=FALSE)
    colori <- c(`filtered out`="khaki", `complete`="powderblue")
    barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
            col = colori, space = 0, main = "Distribution of p-values", ylab="frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
    legend("topright", fill=rev(colori), legend=rev(names(colori)))
    dev.off()
}

sessionInfo()


