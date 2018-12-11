options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(Seurat)
    library(gdata)
    library(optparse)
})

option_list <- list(
    make_option("--data", type="character", help="Input Seurat RDS objectx"),
    make_option("--vars", type="character", default="nUMI", help="Variable(s) used for the regression"),
    make_option("--rds", type="character", help="Output Seurat RDS object")
  )

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

a <- load(args$data)
#change imported object name in seuset if it's not the case
if(!exists("seuset")) mv(a, "seuset")

#Create vector with all variable to remove
list_args = unlist(strsplit(args$vars, ","))
check = list_args %in% colnames(seuset@meta.data)

#check if all variables are in the metadata
if(length(unique(check)) != 1) {
  stop("Some selected variables are not available in your data : ", paste(list_args[check == F], collapse=", "))
}

seuset <- ScaleData(object = seuset, vars.to.regress = list_args)

#Save Seurat object
save(seuset, file = args$rds)

sessionInfo()
