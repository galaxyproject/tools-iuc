# Author: Cico Zhang
# Usage: Rscript bum.R --input p-values.txt --output result.txt --verbose TRUE

# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Avoid crashing Galaxy with an UTF8 error on German LC settings
#loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Import required libraries
suppressPackageStartupMessages({
  library('getopt')
  library('BioNet')
})

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification <- matrix(c(
  'input', 'i', 2, 'character',
  'output', 'o', 2, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options <- getopt(option_specification);

pvals <- read.table(options$input)
bum <- fitBumModel(pvals,plot=FALSE)
mat <- c(bum$lambda, bum$a)
#bumtablename <- paste(options$output,sep="\t")
write.table(x=mat, file=options$output,quote=FALSE, row.names=FALSE, col.names=FALSE)
message ("Parameters have been estimated successfully!")
