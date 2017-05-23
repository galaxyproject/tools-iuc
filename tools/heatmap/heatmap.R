#!/usr/bin/env R
library(argparse)

parser <- ArgumentParser(description='Generate a heatmap from tabular data')

parser$add_argument('--input', dest='input', action="store")
parser$add_argument('--main', dest='main', action="store")
parser$add_argument('--xlab', dest='xlab', action="store")
parser$add_argument('--ylab', dest='ylab', action="store")
parser$add_argument('--var_cols', dest='var_cols', action="append")
parser$add_argument('--scale', dest='scale', action="store")
parser$add_argument('--na_remove', dest='na_remove', action="store")
parser$add_argument('--header', dest='header', action="store")
parser$add_argument('--dendrogram', dest='dendrogram', action="store")
parser$add_argument('--col_min', dest='col_min', action="store")
parser$add_argument('--col_max', dest='col_max', action="store")
parser$add_argument('--output', dest='output', action="store")

args <- parser$parse_args()
options( show.error.messages=F,
         error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) },
         warn=-1 )
inp = read.table( args$input )
x = inp[, as.numeric(unlist(strsplit(args$var_cols, split=",")))]
na_rm_value = FALSE
Colv_value = NA
Rowv_value = NA

if (args$na_remove == 'yes') {
    na_rm_value = TRUE
}

if (args$header == 'yes') {
  colnames(x) = rapply(x[1,], as.character)
  x = x[2:nrow(x),]
}

x = apply(x,2,as.numeric)

if (args$dendrogram == 'row' || args$dendrogram == 'both') {
    Rowv_value = TRUE
}
if (args$dendrogram == 'column' || args$dendrogram == 'both') {
    Colv_value = TRUE
}

## Open output PDF file
pdf( args$output )
rgb_palette = colorRampPalette(c(args$col_min, args$col_max), space="rgb")
heatmap(as.matrix(x), main=args$main, xlab=args$xlab, ylab=args$ylab, scale=args$scale, Rowv=Rowv_value, Colv=Colv_value, na.rm=na_rm_value, col=rgb_palette(10))
## Close the PDF file
devname = dev.off()
