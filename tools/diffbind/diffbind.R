## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library('getopt')
    library('DiffBind')
    library('rjson')
})

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help' , 'h', 0, "logical",
    'outfile' , 'o', 1, "character",
    'plots' , 'p', 2, "character",
    'infile' , 'i', 1, "character",
    'format', 'f', 1, "character",
    'th', 't', 1, "double",
    'bmatrix', 'b', 0, "logical",
    "rdaOpt", "r", 0, "logical"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

parser <- newJSONParser()
parser$addData(opt$infile)
factorList <- parser$getObject()
filenamesIn <- unname(unlist(factorList[[1]][[2]]))
peaks <- filenamesIn[grepl("peaks.bed", filenamesIn)]
bams <- filenamesIn[grepl("bamreads.bam", filenamesIn)]
ctrls <- filenamesIn[grepl("bamcontrol.bam", filenamesIn)]

# get the group and sample id from the peaks filenames
groups <- sapply(strsplit(peaks,"-"), `[`, 1)
samples <- sapply(strsplit(peaks,"-"), `[`, 2)

if ( length(ctrls) != 0 ) {
    sampleTable <- data.frame(SampleID=samples,
                        Condition=groups,
                        bamReads=bams,
                        bamControl=ctrls,
                        Peaks=peaks,
                        stringsAsFactors=FALSE)
} else {
    sampleTable <- data.frame(SampleID=samples,
                        Condition=groups,
                        bamReads=bams,
                        Peaks=peaks,
                        stringsAsFactors=FALSE)
}

cat("sampleTable:\n")
print(sampleTable)

sample = dba(sampleSheet=sampleTable, peakFormat='bed')
sample_count = dba.count(sample)
sample_contrast = dba.contrast(sample_count, categories=DBA_CONDITION, minMembers=2)
sample_analyze = dba.analyze(sample_contrast)
diff_bind = dba.report(sample_analyze, th=opt$th)

# Generate plots
if ( !is.null(opt$plots) ) {
    pdf(opt$plots)
    orvals = dba.plotHeatmap(sample_analyze, contrast=1, correlations=FALSE, cexCol=0.8, th=opt$th)
    dba.plotPCA(sample_analyze, contrast=1, th=opt$th)
    dba.plotMA(sample_analyze, th=opt$th)
    dba.plotVolcano(sample_analyze, th=opt$th)
    dba.plotBox(sample_analyze, th=opt$th)
    dev.off()
}

# Output differential binding sites
resSorted <- diff_bind[order(diff_bind$FDR),]
write.table(as.data.frame(resSorted), file = opt$outfile, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE)

# Output binding affinity scores
if (!is.null(opt$bmatrix)) {
    bmat <- dba.peakset(sample_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
    write.table(as.data.frame(bmat), file="bmatrix.tab", sep="\t", quote=FALSE, row.names=FALSE)
}

# Output RData file
if (!is.null(opt$rdaOpt)) {
    save.image(file = "DiffBind_analysis.RData")
}

sessionInfo()
