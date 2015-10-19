#!/usr/bin/env RScript

# Import required libraries
require('getopt');
library('minfi', quietly=TRUE, warn.conflicts=FALSE,verbose = FALSE)

# Make an option specification, with the following arguments:
# 1. long flag
# 2. short flag
# 3. argument type: 0: No argument, 1: required, 2: optional
# 4. target data type
option_specification = matrix(c(
  'rgset','i',1,'character',
  'preprocess','p',2,'character',
  'cores','c',1,'integer',
), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

# Load the RGset data
if(!is.null(options$rgset)){
	load(options$rgset)
}

# Set phenotype data
pd = pData(options$rgset)
files = gsub(".+/","",pd$filenames)


### Preprocess data with the normalization method chosen

if("${preprocess.preprocess_method}" == "quantile"){
    object = preprocessQuantile(options$rgset)
} else if ("${preprocess.preprocess_method}" == "noob"){
    object = preprocessNoob(options$rgset)
} else if ("${preprocess.preprocess_method}" == "raw"){
    object = preprocessRaw(options$rgset)
} else {
    object = options$rgset
}


M = getM(object)
Beta = getBeta(object)
CN = getCN(object)
chr = seqnames(object)
pos = start(object)


############################################################
# Model Matrix to pass into the bumphunter function
############################################################
pd=pData(object)
T1= as.character(pd$status[2]);T2= as.character(pd$status[1])
keep=pd$status%in%c(T1,T2)
tt=factor(pd$status[keep],c(T1,T2))
design=model.matrix(~tt)

############################################################
# Start bumphunter in a parallel environment
############################################################
# Parallelize over cores on machine
library(doParallel)
registerDoParallel(cores = options$cores)

############################################################
# Bumphunter Run with object processed with Quantile Normalization
############################################################
res=bumphunter(object[,keep],design,B=25,smooth=FALSE,cutoff= 0.3)
bumps= res$tab

############################################################
# Choose DMR's of a certain length threshold.
# This helps reduce the size of DMRs early, and match
# with genes closest to region
############################################################
bumps = bumps[bumps$L>4,]
# gc()
genes=matchGenes(bumps,build="hg19")
result=cbind(bumps,genes)

############################################################
# Save result, which contains DRM's and closest genes
############################################################
save(result,file = "bumps.RData",compress=TRUE)
write.csv(result,file = "bumps.csv")


############################################################
# Garbage collect
############################################################
gc()
