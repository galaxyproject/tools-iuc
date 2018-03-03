#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

                                        # Load libs, common functions, source RaceID,
                                        # Galaxy Params, and read input data (sc)
source(paste(script_dir, "common.R", sep="/"))

                                        # Read input data
message("Count matrix with %.0f cells and %.0f genes", dim(sc@fdata)[1], dim(sc@fdata)[2])
message("Performing clustering using the following parameters: metric=%s, cln=%.0f, do.gap=%.0f clustnr=%.0f, B.gap=%.0f, SE.method=%s, SE.factor=%.2f, bootnr=%.0f, rseed=%.0f", c_metric, c_cln, dogap, c_clustnr, bgap, semethod, sefactor, c_bootnr, c_rseed)

sc <- clustexp(
    sc,
    metric=c_metric, cln=c_cln, do.gap=dogap, clustnr=c_clustnr,
    B.gap=bgap, SE.method=semethod, SE.factor=sefactor,
    bootnr=c_bootnr, rseed=c_rseed
)

message("Plotting images")
plotter("plot_gap", plotgap(sc))
plotter("plot_jaccard", plotjaccard(sc))
plotter("plot_silhouette", plotsilhouette(sc))
plotter("plot_clustheatmap", x <- clustheatmap(sc, final=FALSE, hmethod="single"))


message("Finished plots")

if (generate_final_rdata){
    message("Saving SC object")
    saveRDS(sc, output_rdat)
}





