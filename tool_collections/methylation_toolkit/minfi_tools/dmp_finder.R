#!/usr/bin/Rscript
# Import some required libraries
## Setup R error handling to go to stderr
library('getopt');
### code chunk: Load all required libraries quietly
library('minfi', quietly=TRUE, warn.conflicts=FALSE,verbose = FALSE)

# Make an option specification, with the following arguments:
# 1. long flag
# 2. short flag
# 3. argument type: 0: No argument, 1: required, 2: optional
# 4. target data type
option_specification = matrix(c(
  'rgset','i',1,'character',
  'shrinkVar','s',2,'character'
), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

options$rgset

# Load the RGset data
if(!is.null(options$rgset)){
	load(options$rgset)
}

# Get beta values
beta = getBeta(RGset)

# Set phenotype data
pd = pData(RGset)
phenotype = pd$status

files = gsub(".+/","",pd$filenames)

dmp = dmpFinder(dat=beta,pheno=pd$status,type="categorical",shrinkVar=options$shrinkVar)

write.table(dmp,file="dmpfinder_result.txt",quote=FALSE,row.names=FALSE)