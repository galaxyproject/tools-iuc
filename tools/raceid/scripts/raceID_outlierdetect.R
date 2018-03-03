#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

                                        # Load libs, common functions, source RaceID,
                                        # Galaxy Params, and read input data (sc)
source(paste(script_dir, "common.R", sep="/"))

                                        # Read input data
message("Count matrix with %.0f cells and %.0f genes", dim(sc@fdata)[1], dim(sc@fdata)[2])
message("Detecting outliers using the following parameters: outminc=%.0f, outlg=%.0f, probthr=%.0f outdistquant=%.0f", c_outminc, c_outlg, c_probthr, c_outdistquant)

sc <- findoutliers(
    sc, outminc=c_outminc, outlg=c_outlg, probthr=c_probthr,
    thr=2**-(1:40), outdistquant=c_outdistquant
)

message("Plotting images")
plotter("plot_background", plotbackground(sc))
plotter("plot_sensitivity", plotsensitivity(sc))
plotter("plot_outlierprobs", plotoutlierprobs(sc))
plotter("plot_finalheat", y <- clustheatmap(sc,final=TRUE,hmethod="single"))


message("Finished plots")

message("Generating data tables")
x <- data.frame( CELLID=names(sc@cpart), cluster=sc@cpart )
write.table(
    x[order(x$cluster,decreasing=FALSE),], output_table,
    row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE
)

if (generate_final_rdata){
    message("Saving SC object")
    saveRDS(sc, output_rdat)
}





