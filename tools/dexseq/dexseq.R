## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("DEXSeq")
    library('getopt')
    library('rjson')
})


options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help', 'h', 0, "logical",
    'gtf', 'a', 1, "character",
    'outfile', 'o', 1, "character",
    'reportdir', 'r', 1, "character",
    'rds', 'd', 1, "character",
    'factors', 'f', 1, "character",
    'threads', 'p', 1, "integer",
    'fdr', 'c', 1, "double"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

parser <- newJSONParser()
parser$addData( opt$factors )
factorsList <- parser$getObject()

sampleTable<-data.frame()
countFiles<-c()
factorNames<-c()
primaryFactor<-""
for(factor in factorsList){
    factorName<-factor[[1]]
    factorNames<-append(factorNames, paste(factorName,"exon",sep=":"))
    factorValuesMapList<-factor[[2]]
    c = length(factorValuesMapList)
    for (factorValuesMap in factorValuesMapList){
        for(files in factorValuesMap){
            for(file in files){
                if(primaryFactor == "") {
                    countFiles<-append(countFiles,file)
                }
                sampleTable[basename(file),factorName]<-paste(c,names(factorValuesMap),sep="_")
            }
        }
        c = c-1
    }
    if(primaryFactor == ""){
        primaryFactor <- factorName
    }
}

factorNames<-append(factorNames,"exon")
factorNames<-append(factorNames,"sample")
factorNames<-rev(factorNames)
formulaFullModel <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))
factorNames <- head(factorNames,-1)
formulaReducedModel <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))

sampleTable
formulaFullModel
formulaReducedModel
primaryFactor
countFiles
opt$reportdir
opt$threads
getwd()

dxd = DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design= formulaFullModel, flattenedfile=opt$gtf)

colData(dxd)
dxd <- estimateSizeFactors(dxd)
print("Estimated size factors")
sizeFactors(dxd)
BPPARAM=MulticoreParam(workers=opt$threads)
dxd <- estimateDispersions(dxd, formula=formulaFullModel, BPPARAM=BPPARAM)
print("Estimated dispersions")
dxd <- testForDEU(dxd, reducedModel=formulaReducedModel, fullModel=formulaFullModel, BPPARAM=BPPARAM)
print("tested for DEU")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar=primaryFactor, BPPARAM=BPPARAM)
print("Estimated fold changes")
res <- DEXSeqResults(dxd)
print("Results")
table(res$padj <= opt$fdr)
resSorted <- res[order(res$padj),]
head(resSorted)

export_table <- as.data.frame(resSorted)
last_column <- ncol(export_table)
for(i in 1:nrow(export_table)) {
  export_table[i,last_column] <- paste(export_table[i,last_column][[1]],collapse=", ")
}
write.table(export_table, file = opt$outfile, sep="\t", quote = FALSE, col.names = FALSE)
print("Written Results")

if ( !is.null(opt$rds) ) {
    saveRDS(res, file="DEXSeqResults.rds")
}

if ( !is.null(opt$reportdir) ) {
    DEXSeqHTML(res, fitExpToVar=primaryFactor, path=opt$reportdir, FDR=opt$fdr, color=c("#B7FEA0", "#FF8F43", "#637EE9", "#FF0000", "#F1E7A1", "#C3EEE7","#CEAEFF", "#EDC3C5", "#AAA8AA"))
}
sessionInfo()
