library(devtools)

## Install these
## devtools::install_github(repo = "hhoeflin/hdf5r")
## devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

library(loomR)
library(RaceID)

## Pull an R object

## This R object is NOT an SC object but an object that CONTAINS sc
## Realistically, when moving out of the filter stage, the R object should just be an SC object.
##rdat <- readRDS('out_traject_default.ltree.rdat')

## This is an SCSeq object (specific to RaceID3)
sc <- readRDS('matrix.filter.rdat')

lc <- loomR::create("test3.loom", data=sc@expdata)
#lc <- loomR::connect("test.loom", mode='r+')

## Create a layer that matches the dimensions of the raw_matrix
makeLayer <- function(mat.raw, mat.layer){
    ##
    ## Make a blank matrix of raw size with same names ##
    ##
    mat.new <- matrix(-Inf, nrow=nrow(mat.raw), ncol=ncol(mat.raw))
    rownames(mat.new) <- rownames(mat.raw)
    colnames(mat.new) <- colnames(mat.raw)
    ##
    ## Check that layer hasn't any new genes or cells
    ##
    if (sum(rownames(mat.layer) %in% rownames(mat.new)) != nrow(mat.layer)){
        message("Layer has new features not present in raw_matrix")
    }
    if (sum(colnames(mat.layer) %in% colnames(mat.new)) != ncol(mat.layer)){
        message("Layer has new cells not present in raw_matrix")
    }
    ##
    ## Replace all cell/feature combo sets (slow)
    ##
    for (r in rownames(mat.layer)){
        for (c in colnames(mat.layer)){
            mat.new[r,c] <- mat.layer[r,c]
        }
    }
    return(mat.new)
}

## Add Filter
lc$add.layer(list(filter = makeLayer(sc@expdata, getfdata(sc))))
## Add Normal
lc$add.layer(list(normal = makeLayer(sc@expdata, sc@ndata)))


## generateLoom <- function(sc){
##     loomfile <- loomR::create
##     (
##         filename = filename,
##         data = from@raw.data[, cell.order],
##         cell.attrs = from@meta.data[cell.order, ],
##         layers = list('norm_data' = t(x = from@data[, cell.order])),
##         chunk.dims = chunk.dims,
##         chunk.size = chunk.size,
##         overwrite = overwrite,
##         display.progress = display.progress
##     )
## }
