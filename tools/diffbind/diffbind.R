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
    'infile' , 'i', 1, "character",
    'outfile' , 'o', 1, "character",
    'scorecol', 'n', 1, "integer",
    'lowerbetter', 'l', 1, "logical",
    'summits', 's', 1, "integer",
    'th', 't', 1, "double",
    'format', 'f', 1, "character",
    'plots' , 'p', 2, "character",
    'bmatrix', 'b', 0, "logical",
    "rdaOpt", "r", 0, "logical",
    'infoOpt' , 'a', 0, "logical",
    'verbose', 'v', 2, "integer",
    'help' , 'h', 0, "logical"
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
                        Tissue=samples, # using "Tissue" column to display ids as labels in PCA plot
                        stringsAsFactors=FALSE)
} else {

    sampleTable <- data.frame(SampleID=samples,
                        Replicate=samples,
                        Condition=groups,
                        bamReads=bams,
                        Peaks=peaks,
                        Tissue=samples,
                        stringsAsFactors=FALSE)
}

sample = dba(sampleSheet=sampleTable, peakFormat='bed', scoreCol=opt$scorecol, bLowerScoreBetter=opt$lowerbetter)

if ( !is.null(opt$summits) ) {
    sample_count = dba.count(sample, summits=opt$summits)
} else {
    sample_count = dba.count(sample)
}

sample_contrast = dba.contrast(sample_count, categories=DBA_CONDITION, minMembers=2)
sample_analyze = dba.analyze(sample_contrast)
diff_bind = dba.report(sample_analyze, th=opt$th)

# Generate plots
if ( !is.null(opt$plots) ) {
    pdf(opt$plots)
    orvals = dba.plotHeatmap(sample_analyze, contrast=1, correlations=FALSE, cexCol=0.8, th=opt$th)
    dba.plotPCA(sample_analyze, contrast=1, th=opt$th, label=DBA_TISSUE, labelSize=0.3)
    dba.plotMA(sample_analyze, th=opt$th)
    dba.plotVolcano(sample_analyze, th=opt$th)
    dba.plotBox(sample_analyze, th=opt$th)
    dev.off()
}

# Output differential binding sites
resSorted <- diff_bind[order(diff_bind$FDR),]
# Convert from GRanges (1-based) to 0-based format (adapted from https://www.biostars.org/p/89341/)
if (opt$format == "bed") {
    resSorted  <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)))
} else if (opt$format == "interval") {
     # Output as interval
    df <- as.data.frame(resSorted)
    extrainfo <- NULL
    for (i in 1:nrow(df)) {
        extrainfo[i] <- paste0(c(df$width[i], df[i, 6:ncol(df)]), collapse="|")
    }
    resSorted  <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)),
        Comment=extrainfo)
} else {
    # Output as 0-based tabular
    resSorted <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)),
        mcols(resSorted))
}
write.table(resSorted, file = opt$outfile, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE)

# Output binding affinity scores
if (!is.null(opt$bmatrix)) {
    bmat <- dba.peakset(sample_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
    # Output as 0-based tabular
    bmat <- data.frame(Chrom=bmat[, 1],
        Start=bmat[, 2] - 1,
        End=bmat[, 3],
        bmat[, 4:ncol(bmat)])
    write.table(bmat, file="bmatrix.tab", sep="\t", quote=FALSE, row.names=FALSE)
}

# Output RData file
if (!is.null(opt$rdaOpt)) {
    save.image(file = "DiffBind_analysis.RData")
}

# Output analysis info
if (!is.null(opt$infoOpt)) {
    info <- "DiffBind_analysis_info.txt"
    cat("dba.count Info\n\n", file=info, append = TRUE)
    capture.output(sample, file=info, append=TRUE)
    cat("\ndba.analyze Info\n\n", file=info, append = TRUE)
    capture.output(sample_analyze, file=info, append=TRUE)
    cat("\nSessionInfo\n\n", file=info, append = TRUE)
    capture.output(sessionInfo(), file=info, append=TRUE)
}