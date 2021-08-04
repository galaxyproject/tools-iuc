
suppressWarnings(suppressPackageStartupMessages(library(xbioc)))
suppressWarnings(suppressPackageStartupMessages(library(MuSiC)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))

args = commandArgs(trailingOnly=TRUE)
source(args[1])

if (dataset == "bulk"){
    rds = bulk_eset
} else {
    rds = scrna_eset
}

printout <- function(text){
    if (typeof(text) %in% c("list","vector")){
        write.table(text, file = outfile_tab, quote=F, sep="\t", col.names=NA)
    } else {
        ## text
        capture.output(text, file = outfile_tab)
    }    
}

if (inspector %in% c("pData", "fData", "dims", "experimentData", "signature", "annotation", "abstract")){
    op = get(inspector)
    tab = op(rds)
    printout(tab)    
} else {
    stop(paste0("No such option:", inspector))
}