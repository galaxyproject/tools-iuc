##!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

                                        # Load libs, common functions, source RaceID,
                                        # Galaxy Params, and read input data (sc)
source(paste(script_dir, "common.R", sep="/"))

                                        # Read input data
message("Count matrix with %.0f cells and %.0f genes", dim(sc@fdata)[1], dim(sc@fdata)[2])

sc <- comptsne(sc,rseed = c_rseed)

message("Plotting initial and final tSNEs")

png("plot_initial.png")
a <- plottsne(sc,final = F)
dev.off()
png("plot_final.png")
a <- plottsne(sc,final = T)
dev.off()

if (gene_sets != "" ){
####all_sets <- strsplit(gene_sets, "+?\\s*__split__\\s*,?")
    all_sets <- strsplit(gene_sets, '\\s*\\+?\\s*_split_\\s*,?')
    print(all_sets)

    for (given in all_sets){
        given <- trimws(given)
        message("Plotting %s", given)
        g <- c(unlist(strsplit(given,'\\s*\\+\\s*')))
        png(paste("plot", given, sep="_"))
####print(plotexptsne(sc,g))
        plotexptsne(sc,g, n=given, logsc=T)
        dev.off()
    }
}


if (regex_val != ""){
    message("using subcell groups")
    png("plot_labels.png")
    plotlabelstsne(sc, labels = sub(regex_val, "", names(sc@ndata)))
    dev.off()
    png("plot_symbols.png")
    plotsymbolstsne(sc, types = sub(regex_val, "", names(sc@ndata)))
    dev.off()
} else {
    message("using all cell ids")
    png("plot_labels.png")
    plotlabelstsne(sc, labels = names(sc@ndata) )
    dev.off()
    png("plot_symbols.png")
    plotsymbolstsne(sc, types = names(sc@ndata) )
    dev.off()
}


message("Saving SC object")
saveRDS(sc, output_rdat)
