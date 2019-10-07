options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(GWASTools)
    library(optparse)
})
option_list <- list(
    make_option(c("-f", "--file"), type="character", help="RInput GWAS file"),
    make_option("--pval", type="integer", help="Pvalue column"),
    make_option("--chromosome", type="integer", help="Chromosome column"),
    make_option("--ymin", type="double", help="Min y value"),
    make_option("--ymax", type="double", help="Max y value"),
    make_option("--trunc", help="Show truncation lines", action="store_true"),
    make_option("--sig", type="double", help="Significance level for lines"),
    make_option("--thin", type="double", help="Thinning value", action="store_true", dest="thin"),
    make_option("--ppb", type="integer", help="Points per bin, if thinning value is specified", action="store_true"))
args <- parse_args(OptionParser(option_list=option_list))
file <- args$file
data <-  read.table(args$file, header=TRUE)
pval <-  data[,args$pval]
chromosome <-  data[,args$chromosome]
if(!is.null(args$ymin) & !is.null(args$ymax)){
    ylimit <-  c(args$ymin,args$ymax)
}else if(xor(!is.null(args$ymin), !is.null(args$ymax))){
    print("If specifying range, both ymin and ymax must be set")
    ylimit <- NULL
}else{
    ylimit <- NULL
}
if(is.null(args$trunc)){
    trunc <- FALSE
}else{
    trunc <- TRUE
}
if(!is.null(args$sig)){
    sig <-  args$sig
}else{
    sig <- NULL
}
if(!is.null(args$thin)){
    thin = args$thin
    if(thin == 0){
        thin = NULL 
    }
}else{
    thin <- FALSE
}
if(!is.null(thin) & !is.null(args$ppb)){
    ppb = args$ppb
}
pdf("manhattan.pdf")
if(isFALSE(thin)){
    manhattanPlot(pval, chromosome, ylim=ylimit, trunc.lines=trunc, signif=sig)
}else{
    manhattanPlot(pval, chromosome, ylim=ylimit, trunc.lines=trunc, signif=sig, thinThreshold = thin, pointsPerBin = ppb)   
}
invisible(dev.off())
