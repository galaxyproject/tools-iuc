#!/usr/bin/env R
library(argparse)
library(RColorBrewer)

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
parser$add_argument('--palette', dest='palette', action="store")
parser$add_argument('--shades', dest='shades', action="store")
parser$add_argument('--select_color', dest='select_color', action="store_true")

args <- parser$parse_args()

options( show.error.messages=F,
         error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) },
         warn=-1 )

header = FALSE
if (args$header == 'yes') {
  header = TRUE
}

inp = as.matrix(read.table(args$input, header=header, row.names=1, as.is=TRUE, check.names=FALSE))

columns = lapply(as.numeric(unlist(strsplit(args$var_cols, split=","))), function(x) x - 1)

x = apply(inp[1:nrow(inp),1:ncol(inp)], c(1,2), as.numeric)[,unlist(columns)]

na_rm_value = FALSE
Colv_value = NA
Rowv_value = NA

if (args$na_remove == 'yes') {
    na_rm_value = TRUE
}

if (args$dendrogram == 'row' || args$dendrogram == 'both') {
    Rowv_value = TRUE
}

if (args$dendrogram == 'column' || args$dendrogram == 'both') {
    Colv_value = TRUE
}

if (args$select_color) {
    rgb_palette = colorRampPalette(c(args$col_min, args$col_max), space="rgb")
    colors = rgb_palette(as.numeric(args$shades))
} else {
    colors = brewer.pal(as.numeric(args$shades), args$palette)
}

## Open output PDF file
pdf( args$output )
heatmap(x, main=args$main, xlab=args$xlab, ylab=args$ylab, scale=args$scale, Rowv=Rowv_value, Colv=Colv_value, na.rm=na_rm_value, col=colors)
## Close the PDF file
devname = dev.off()
