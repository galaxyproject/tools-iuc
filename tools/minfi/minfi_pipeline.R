# TODO: Remove print statements
# TODO: Make sure to Garbage collect

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
    "quiet", "q", 0, "logicall",
    "help", "h", 0, "logical",
    "preprocess","p",2,"character",
    "cores","c",1,"integer")
    ,byrow=TRUE, ncol=4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

# enforce the following required arguments

if (is.null(opt$preprocess)) {
    cat("'preprocess' is required\n")
    q(status=1)
}

verbose <- if (is.null(opt$quiet)) {
    TRUE
} else {
    FALSE
}

# Load required libraries

suppressPackageStartupMessages({
    library("minfi")
    library("FlowSorted.Blood.450k")
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    library("doParallel")
})


## Parse cheetah code and make dataframe for creating tmp dir
minfi_config_file = paste0("minfi_temp","/minfi_config.txt")
minfi_config = read.table(minfi_config_file)
colnames(minfi_config) = c("status","green","red","name")

minfi_config

## Make the tmpdir for symlinking data
base_dir = paste0("minfi_temp","/base")
system(paste0("mkdir ",base_dir))


## Make symlinks of files
for (i in 1:nrow(minfi_config)){
    stopifnot(nrow(minfi_config) == nrow(minfi_config["name"]))

    ## Make green idat file symlinks
    file_green = paste0(base_dir,"/",as.character(minfi_config[i,"name"]),"_Grn.idat")
    cmd_green = paste("ln -s",as.character(minfi_config[i,"green"]),file_green,sep=" ")
    cat("Reading file ",i,"GREEN Channel ", file_green)
    system(cmd_green)

    ## Make red idat file symlinks
    file_red = paste0(base_dir,"/",as.character(minfi_config[i,"name"]),"_Red.idat")
    cmd_red = paste("ln -s",as.character(minfi_config[i,"red"]),file_red,sep=" ")
    cat("Reading file ",i,"RED Channel ", file_red)
    system(cmd_red)
}

## Make dataframe with Basenames
Basename = paste0(base_dir,"/",unique(substr(list.files(base_dir),1,17)))
status = minfi_config[match(gsub(".+/","",Basename), minfi_config$name),"status"]
targets = data.frame(Basename=Basename,status=status)

targets

## Read 450k files
RGset = read.450k.exp(targets=targets,verbose=T)
RGset


## Preprocess data with the normalization method chosen
if(opt$preprocess == "quantile"){
    normalized_RGset = preprocessQuantile(RGset)
    print("Data has been normalized using Quantile Normalization")
} else if (opt$preprocess == "noob"){
    normalized_RGset = preprocessNoob(RGset)
    print("Data has been normalized using Noob Normalization")
} else if (opt$preprocess == "raw"){
    normalized_RGset = preprocessRaw(RGset)
    print("Data has been normalized using Raw Normalization")
} else if (opt$preprocess == "illumina"){
    normalized_RGset = preprocessIllumina(RGset,bg.correct = TRUE, normalize = c("controls", "no"),reference = 1)
    print("Data has been normalized using Illumina Normalization")
} else if (opt$preprocess == "preprocessFunnorm"){
    normalized_RGset = preprocessFunnorm(RGset)
    print("Data has been normalized using Functional Normalization")
} else {
    normalized_RGset = RGset
    print("No Normalization applied")
}


## Get beta values from Proprocessed data
beta = getBeta(normalized_RGset)

## Set phenotype data
pd = pData(RGset)


## QC REPORT
files = gsub(".+/","",pd$filenames)
## Produce PDF file
if (!is.null(RGset)) {
    # Make PDF of QC report
    minfi::qcReport(rgSet=RGset,sampNames=files,sampGroups=pd$status,pdf="qc_report.pdf")
}


## MDS Plot
## Set phenotype data
files = gsub(".+/","",pd$filenames)
#numPositions=as.integer("${numPositions}")

## Produce PDF file
if (!is.null(RGset)) {
    ## Make PDF of MDS plot
    pdf("mds_plot.pdf")
    minfi::mdsPlot(dat=RGset,sampNames=files,sampGroups=pd$status,main="Beta MDS",numPositions=1000,pch=19)
    dev.off()
}

## TODO: Fix Estimate cell counts!
# Estimate Cell counts
#result = minfi::estimateCellCounts(dat=RGset,meanPlot="${meanplot}")
#write.table(result,file="estimated_cell_counts.txt",quote=FALSE,row.names=FALSE)

print("Plot making finished")

## DMP finder
dmp = dmpFinder(dat=beta,pheno=pd$status,type="categorical")
#,shrinkVar="${shrinkVar}")
write.table(dmp,file="dmps.csv",quote=FALSE,row.names=TRUE)
print("DMP Finder successful")


#M = getM(normalized_RGset)
#Beta = getBeta(normalized_RGset)
#CN = getCN(normalized_RGset)
#chr = seqnames(normalized_RGset)
#pos = start(normalized_RGset)

#print("retrieved Beta and meth values")

# Model Matrix to pass into the bumphunter function
pd=pData(normalized_RGset)
T1= levels(pd$status)[2]
T2= levels(pd$status)[1]

# Introduce error if levels are different
stopifnot(T1!=T2)

keep=pd$status%in%c(T1,T2)
tt=factor(pd$status[keep],c(T1,T2))
design=model.matrix(~tt)

design

# Start bumphunter in a parallel environment
# Parallelize over cores on machine
registerDoParallel(cores = opt$cores)

## Bumphunter Run with normalized_RGset processed with Quantile Normalization
res=bumphunter(normalized_RGset[,keep],design,B=25,smooth=FALSE,cutoff= 0.3,type="Beta")
bumps= res$tab

print("DMRs found with bumphunter")
head(bumps)

## Choose DMR's of a certain length threshold.
## This helps reduce the size of DMRs early, and match
## with genes closest to region
print("Excluding DMR's with length smaller than 4")
bumps = bumps[bumps$L>4,]

# Save result, which contains DMR's and closest genes
write.table(bumps,file = "dmrs.csv",quote=FALSE)

# Garbage collect
gc()


## TODO: FIX BLOCK FINDER
# Block finder
#library(sva)
#pheno <- pData(GRset)
#mod <- model.matrix(~as.factor(status), data=pheno)
#mod0 <- model.matrix(~1, data=pheno)
#sva.results <- sva(mval, mod, mod0)
