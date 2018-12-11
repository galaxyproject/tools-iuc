options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(optparse)
    library(gdata)
})

option_list <- list(
    make_option("--data", type="character", help="Input Seurat RDS objectx"),
    make_option("--scale", type="double", default=NULL, help="Scale factor used for normalization"),
    make_option("--rds", type="character", help="Output Seurat RDS object")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

a <- load(args$data)
#change imported object name in seuset if it's not the case
if(!exists("seuset")) mv(a, "seuset")

if(is.null(args$scale)){
    scale_factor = median(seuset@meta.data$nUMI)
} else {
    scale_factor = args$scale
}

seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", 
    scale.factor = scale_factor)

#Save Seurat object
save(seuset, file = args$rds)

sessionInfo()
