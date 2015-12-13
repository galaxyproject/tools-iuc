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
    "preprocess","p",1,"character",
    "cores","c",1,"integer",
    "numPositions","n",2,"integer",
    "shrinkVar","s",2,"logical",
    "b_permutations","b",2,"integer",
    "smooth","m",2,"logical",
    "cutoff","t",2,"double",
    "l_value","l",2,"integer")
    ,byrow=TRUE, ncol=4)
opt <- getopt(spec)

# If help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}


## Set verbose mode
verbose = if(is.null(opt$quiet)){TRUE}else{FALSE}
if(verbose){
    cat("Verbose mode is ON\n\n")
}

# Enforce the following required arguments
if (is.null(opt$preprocess)) {
    cat("'--preprocess' is required\n")
    q(status=1)
}
cat("preprocess = ",opt$preprocess,"\n")
cat("cores = ", opt$cores, "\n")
cat("b_permutations = ",opt$b_permutations,"\n")
cat("smooth = ",opt$smooth,"\n")
cat("cutoff = ",opt$cutoff,"\n")
cat("l_value = ",opt$l_value,"\n")
cat("numPositions = ",opt$numPositions,"\n")
cat("shrinkVar = ",opt$shrinkVar,"\n")


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

if ( verbose ) {
    cat("Minfi targets file:\n\n ")
    print(targets)
}

## Read 450k files
RGset = read.450k.exp(targets=targets,verbose=T)

if (verbose){
    cat("RGset has been read: \n\n")
    print(RGset)
}

pd = pData(RGset)

## NOTE: QC report is for samples before normalization
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
    minfi::mdsPlot(dat=RGset,sampNames=files,sampGroups=pd$status,main="Beta MDS",numPositions = opt$numPositions,pch=19)
    dev.off()
}

if(verbose){
    cat("Made plot of QC report and MDS plot\n\n")
}


## Preprocess data with the normalization method chosen
if(opt$preprocess == "quantile"){
    normalized_RGset = preprocessQuantile(RGset)
    if (verbose){cat("Preprocessed using Quantile normalization")};
} else if (opt$preprocess == "funnorm"){
    normalized_RGset = preprocessFunnorm(RGset)
    if(verbose){print("Preprocessed using Functional normalization")}
} else if (opt$preprocess == "noob"){
    normalized_RGset = preprocessNoob(RGset)
    if (verbose){cat("Preprocessed using Noob normalization")};
} else if (opt$preprocess == "illumina"){
    normalized_RGset = preprocessIllumina(RGset,bg.correct = TRUE, normalize = c("controls", "no"),reference = 1)
    if(verbose){print("Preprocessed using Illumina normalization")}
} else if (opt$preprocess == "swan"){ 
    normalized_RGset = preprocessSWAN(RGset)
    if(verbose){print("Preprocessed using SWAN normalization")}
}else {
    normalized_RGset = RGset
    if(verbose){print("Preprocessed using No normalization")}
}



## Get beta values from Proprocessed data
beta = getBeta(normalized_RGset)
## Set phenotype data
pd = pData(normalized_RGset)


## DMP finder
dmp = dmpFinder(dat=beta,pheno=pd$status,type="categorical",shrinkVar=opt$shrinkVar)
write.csv(dmp,file="dmps.csv",quote=FALSE,row.names=TRUE)
if(verbose){
    cat("\n","DMP Finder successful \n")
}


# Model Matrix to pass into the bumphunter function
T1= levels(pd$status)[2]
T2= levels(pd$status)[1]

# Introduce error if levels are different
stopifnot(T1!=T2)

keep=pd$status%in%c(T1,T2)
tt=factor(pd$status[keep],c(T1,T2))
design=model.matrix(~tt)

if(verbose){
    cat("Model matrix is: \n")
    design
}

# Start bumphunter in a parallel environment
# Parallelize over cores on machine
registerDoParallel(cores = opt$cores)

## Bumphunter Run with normalized_RGset processed with Quantile Normalization
if (is(normalized_RGset,"GenomicRatioSet")) {
    res=bumphunter(normalized_RGset[,keep],design,B=opt$b_permutations,smooth=opt$smooth,cutoff= opt$cutoff,type="Beta")
    bumps= res$tab
} else if(is(normalized_RGset,"MethylSet")) {
    # convert MethylSet (norm.swan) into GenomicRatioSet through Ratioset
    normalized_RGset = ratioConvert(normalized_RGset, what = "both", keepCN = TRUE)
    normalized_RGset = mapToGenome(normalized_RGset)
    res = bumphunter(normalized_RGset[,keep],design=design,B=opt$b_permutations,smooth=opt$smooth,cutoff=opt$cutoff)
    bumps = res$tab
} else {
    # This else case is never supposed to run,
    # it will run only if the normalized_RGset object
    # did not have the expected class type   
    cat("Bumphunter did not run properly","\n")
    stopifnot(!is.null(bumps))
}

if(verbose){
    cat("Bumphunter result", "\n")
    head(bumps)
}


## Choose DMR's of a certain length threshold.
## This helps reduce the size of DMRs early, and match
## with genes closest to region
bumps = bumps[bumps$L > opt$l_value,]
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tab=matchGenes(bumps,genes)
annotated_dmrs=cbind(bumps,tab)

if(verbose){
    cat("Match with annotation\n")
    head(annotated_dmrs)
}

# Save result, which contains DMR's and closest genes
write.csv(annotated_dmrs,file = "dmrs.csv",quote=FALSE,row.names=FALSE)

# Garbage collect
gc()


## TODO: FIX BLOCK FINDER
# Block finder
#library(sva)
#pheno <- pData(GRset)
#mod <- model.matrix(~as.factor(status), data=pheno)
#mod0 <- model.matrix(~1, data=pheno)
#sva.results <- sva(mval, mod, mod0)
