#!/usr/bin/env Rscript

suppressMessages(library(optparse, warn.conflicts = FALSE))
opt_list=list(
make_option(c("-d", "--directory"), type="character", default=NULL, help="directory containing the samples", metavar="character"),
make_option(c("-p", "--phendat"), type="character", default=NULL, help="phenotype data(must be a .csv file)", metavar="character"),
make_option(c("-t","--outputtranscript"), type="character", default="output_transcript.csv", help="output_transcript.csv: contains the transcripts of the experiments", metavar="character"),
make_option(c("-g","--outputgenes"), type="character", default="output_genes.csv", help="output_genes.csv: contains the genes of the experiments", metavar="character"),
make_option(c("-e","--texpression"), type="double", default="0.5", help="transcripts expression filter", metavar="character"),
make_option(c("--bgout"), type="character", default="", help="save the ballgown object created in the process", metavar="character"),
make_option(c("-f","--format"), type="character", default="tsv", help="Create csv or tsv files as output", metavar="character"),
make_option(c("-T","--tsvoutputtranscript"), type="character", default="output_transcript.tsv", help="output_transcript.tsv: contains the transcripts of the experiments", metavar="character"),
make_option(c("-G","--tsvoutputgenes"), type="character", default="output_genes.tsv", help="output_genes.tsv: contains the genes of the experiments", metavar="character")
)
opt_parser=OptionParser(option_list=opt_list)
opt=parse_args(opt_parser)

# Loading required libraries. suppressMessages() remove all noisy attachement messages
# ----------------------------------------------------------------------------------------

suppressMessages(library(ballgown, warn.conflicts = FALSE))
suppressMessages(library(genefilter, warn.conflicts = FALSE))
suppressMessages(library(dplyr, warn.conflicts = FALSE))

# Setup for the tool with some bases variables.
# ----------------------------------------------------------------------------------------


filtstr = opt$texpression
pdat = 2
phendata = read.csv(opt$phendat)


# Checking if the pdata file has the right samples names.
# ----------------------------------------------------------------------------------------

if (all(phendata$ids == list.files(opt$directory)) != TRUE)
{
  stop("Your phenotype data table does not match the samples names. ")
}

# Creation of the ballgown object based on data
# ----------------------------------------------------------------------------------------
bgi = ballgown(dataDir= opt$directory , samplePattern="", pData = phendata, verbose = FALSE)

# Filter the genes with an expression superior to the input filter
# ----------------------------------------------------------------------------------------
bgi_filt= subset(bgi, paste("rowVars(texpr(bgi)) >",filtstr), genomesubset = TRUE)

# Creating the variables containing the transcripts and the genes and sorting them through the arrange() command.
# Checking if there's one or more adjust variables in the phenotype data file
# ----------------------------------------------------------------------------------------

if (ncol(pData(bgi))<=3) {
  results_transcripts=stattest(bgi_filt, feature = "transcript", covariate = colnames(pData(bgi))[pdat], adjustvars = colnames(pData(bgi)[pdat+1]), getFC = TRUE, meas = "FPKM")
  results_genes=stattest(bgi_filt, feature = "gene", covariate = colnames(pData(bgi))[pdat], adjustvars = colnames(pData(bgi)[pdat+1]), getFC = TRUE, meas = "FPKM")
} else {
  results_transcripts=stattest(bgi_filt, feature = "transcript", covariate = colnames(pData(bgi))[pdat], adjustvars = c(colnames(pData(bgi)[pdat+1:ncol(pData(bgi))])), getFC = TRUE, meas = "FPKM")
  results_genes=stattest(bgi_filt, feature = "gene", covariate = colnames(pData(bgi))[pdat], adjustvars = c(colnames(pData(bgi)[pdat+1:ncol(pData(bgi))])), getFC = TRUE, meas = "FPKM")
}

results_transcripts = data.frame(geneNames=ballgown::geneNames(bgi_filt), geneIDs=ballgown::geneIDs(bgi_filt), results_transcripts)
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)

# Main output of the wrapper, two .csv files containing the genes and transcripts with their qvalue and pvalue
#This part also output the data of the ballgown object created in the process and save it in a R data file
# ----------------------------------------------------------------------------------------
if (opt$format == "tsv"){
  write.table(results_transcripts, file=opt$tsvoutputtranscript, quote=FALSE, sep='\t', col.names = NA)
  write.table(results_genes, file=opt$tsvoutputgenes, quote=FALSE, sep='\t', col.names = NA)
} else {
  write.csv(results_transcripts, opt$outputtranscript, row.names = FALSE)
  write.csv(results_genes, opt$outputgenes, row.names = FALSE)
}

if (opt$bgout != ""){
  save(bgi, file=opt$bgout)
}
