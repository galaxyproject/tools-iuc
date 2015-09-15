###################################################
# code chunk number 1: Record starting time
###################################################
timeStart <- as.character(Sys.time())

### code chunk number 2: Load all required libraries
msg.out1 <- capture.output( suppressMessages( library(IlluminaHumanMethylation450kanno.ilmn12.hg19,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)))
msg.out2 <- capture.output( suppressMessages( library(minfi,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)))
require("getopt",quietly=TRUE)

### code chunk number 3: Update package version
if (packageVersion("minfi") < "1.15.6") {
  stop("Please update 'minfi' to version >= 1.15.7  to run this tool")
}

# Make an option specification, with the following arguments:
# 1. input is a gdac downloaded TCGA matrix
# 2. output is an optional

option_specification = matrix(c(
  'input', 'i', 1, 'character',
  'output','o', 2, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# Read the TCGA data
dat = readTCGA(options$input, sep = "\t", keyName = "Composite Element REF", Betaname = "Beta_value", pData = NULL, array = "IlluminaHumanMethylation450k",showProgress=FALSE)

#Get beta values
beta = getBeta(dat)

# Store beta values in a
write.table(beta, file = options$output ,row.names = F,quote = F, sep="\t")
