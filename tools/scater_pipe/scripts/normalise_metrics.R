#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

## Load libs, common functions, source Galaxy config
source(paste(script_dir, "common.R", sep="/"))

require(scater)

## Already have list 
tab2html <- function(tabl){
    str = "<table>"
    str = paste(str, "<tr>", sep="\n")
    str = paste(str, "<th>", "Method", "</th>", sep="")
    for (header in colnames(tabl)){
        str = paste(str, "<th>", header, "</th>", sep="")
    }
    str = paste(str, "</tr>", sep="")
    
    # rows    
    r = 0
    while (r < dim(tabl)[1]){
        r <- r + 1
        str <- paste(str,
            "<tr><td><b>", rownames(tabl)[r], 
                     "</b></td><td>",
            paste(signif(tab[r,],6), collapse="</td><td>"),
            "</td></tr>",
            sep = ""
        )
    }
    str <- paste(str, "</table>", sep="\n")
    return(str)
}

## We have matrices, which is a map of matrix names and objects
metric <- logcounts

exprslog <- c()
tableau <- c()
exprsmat <- list()

for (name in names(matrices)){
    exprsmat[[name]] <- metric(matrices[[name]])
    mat <- exprsmat[[name]]
    first_line <- head(mat,1)
    tableau <- rbind(tableau, first_line)
    exprslog <- c(exprslog, T)
}

rownames(tableau) <- names(matrices)
tab <- tableau[,1:6]

                                        #plotter("OUT_plotrle" , function(){
svg("OUT_plotrle.svg", width=10, height=10)
plotRLE(
    matrices[[1]],
    exprs_mats = exprsmat,
    exprs_logged = exprslog,
    colour_by = colourby,
    order_by_color = FALSE
)
dev.off()

write(tab2html(tab), "table.html")
