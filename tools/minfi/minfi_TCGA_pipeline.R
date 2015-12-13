# setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    'quiet', 'q', 2, "logical",
    'help' , 'h', 0, "logical",
    'tarfile','f',1,"character",
    "cores","c",1,"integer",
    "b_permutations","b",2,"integer",
    "smooth","m",2,"logical",
    "l_value","l",2,"integer")
    ,byrow=TRUE, ncol=4)
opt <- getopt(spec)

## If help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}


## Set verbose mode
verbose = if(is.null(opt$quiet)){TRUE}else{FALSE}
if(verbose){
    cat("Verbose mode is ON\n\n")
}

## Load required libraries
suppressPackageStartupMessages({
    library("minfi")
    library("FlowSorted.Blood.450k")
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    library("doParallel")
    library("tools")
})


config_file = "tcga_temp/config.txt"
conf = read.csv(config_file,stringsAsFactors=FALSE,header=F)
tarfile_name = gsub(".+/","",conf$V2)
dataset_path = conf$V1

cmd = paste("ln -s",dataset_path,tarfile_name,sep=" ")
cat("Command : ", cmd,"\n")
system(cmd)

tarfile = tarfile_name 
cat ("tarfile name: ",tarfile," file ext: ",file_ext(tarfile))
## UNtar files in R first
if (file_ext(tarfile) == "tar"){
    cat("Entering IF statment")
    tar_contents = untar(tarfile,list=TRUE)
    cat("regex failing here")
    f = as.character(tar_contents[grep(".data.txt",fixed=TRUE,x=tar_contents)])
    if (!is.null(f)){
        cat("Untar being attempted")
        untar(tarfile,  exdir = ".",files=f)
        cat("Untar succcess")
    }
}
## Move file from sub directory to main directory
from = list.files(pattern=".data.txt",recursive=TRUE)
to = gsub(".+/","",from)
rename_success = file.rename(from=from, to=to)    


# This should pass only if steps have been successful 
stopifnot(rename_success)
if (rename_success){
    input_file = to 
}

### Read the TCGA data
GRset = readTCGA(input_file, sep = "\t", keyName = "Composite Element REF", Betaname = "Beta_value", pData = NULL, array = "IlluminaHumanMethylation450k",showProgress=TRUE)

### Get beta values
beta = getBeta(GRset)
pd = pData(GRset)
CN = getCN(GRset)
chr = seqnames(GRset)
pos = start(GRset)
chr = as.character(chr)


# Assign phenotype information
## Based on TCGA sample naming, TCGA-2E-A9G8-01A-11D-A409-05, char 14,15 represent
## phenotypic status of sample, 01 = cancer, 11=normal 
pd$status = ifelse(test= (substr(rownames(pd),14,15) == "01"),yes="cancer",no="normal")

### DMP finder
dmp = dmpFinder(dat=beta,pheno=pd$status,type="categorical")
write.csv(dmp,file="dmps.csv",quote=FALSE,row.names=TRUE)
if(verbose){
    cat("DMP Finder successful \n")
}

 
## Make design matrix for bumphunting 
#Model Matrix
T1="normal";T2="cancer"
## Introduce error if levels are different
stopifnot(T1!=T2)
keep=pd$status%in%c(T1,T2)
tt=factor(pd$status[keep],c(T1,T2))
design=model.matrix(~tt)

 
## Start bumphunter in a parallel environment
## Parallelize over cores on machine
library(doParallel)
registerDoParallel(cores = opt$cores)

# Bumphunter Run with object processed with default Quantile Normalization
# provided along with TCGA data
dmrs = bumphunter(beta[,keep],chr=chr,pos=pos,design=design,B=opt$b_permutations,smooth=opt$smooth,pickCutoff =TRUE)
bumps = dmrs$tab

if(verbose){
    cat("Bumphunter result", "\n")
    head(bumps)
}


### Choose DMR's of a certain length threshold.
### This helps reduce the size of DMRs early, and match
### with genes closest to region
bumps = bumps[bumps$L > opt$l_value,]
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tab=matchGenes(bumps,genes)
annotated_dmrs=cbind(bumps,tab)

if(verbose){
    cat("Match with annotation\n")
    head(annotated_dmrs)
}

## Save result, which contains DMR's and closest genes
write.csv(annotated_dmrs,file = "dmrs.csv",quote=FALSE,row.names=FALSE)

## Garbage collect
gc()


