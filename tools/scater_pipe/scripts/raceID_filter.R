#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

# Load libs, common functions, source RaceID and Galaxy Params
source(paste(script_dir, "common.R", sep="/"))

# Read input data
x <- read.csv(count_matrix, sep="\t", header=TRUE)
rownames(x) <- x[,1]

message("Count matrix with %.0f cells and %.0f genes", dim(x)[1], dim(x)[2])

# Control gene filtering
# (if blank do nothing)
if (!( (!exists("control_genes_filter")) || control_genes_filter == "" || control_genes_filter == "None")){
    c_genes <- unlist(strsplit(control_genes_filter, "\\s*,\\s"))
    for (cg in c_genes){
        x <- x[grep(cg, rownames(x), invert=T), -1]
        message("Filtering against %s yielded %.0f cells and %.0f genes", cg, dim(x)[1], dim(x)[2])
    }
} else {
  x <- x[grep("INVALID", rownames(x), invert=T), -1]
}

if (!filtering){
   # No filtering, just return an SCseq object
   sc <- SCseq(x)
   saveRDS(sc, output_rdat)
   quit(status=0)
}


sc <- SCseq(x)

if (c_maxexpr == 0){
   c_maxexpr = Inf
}

# Perform actual filtering beyond this point
message("Applying filtering parameters: mintotal = %.0f, minexpr = %.0f, minnumber= %.0f, maxexpr= %.0f, downsample= %.0f, dsn= %.0f, rseed=%.0f", c_mintotal, c_minexpr, c_minnumber, c_maxexpr, c_downsample, c_dsn, c_rseed)

sc <- filterdata(sc,
   mintotal=c_mintotal,
   minexpr=c_minexpr,
   minnumber=c_minnumber,
   maxexpr=c_maxexpr, downsample = c_downsample, dsn=c_dsn, rseed=c_rseed)
   
message("Output matrix yielded %.0f cells and %.0f genes", dim(sc@fdata)[1], dim(sc@fdata)[2])

# Output table
write.table(sc@fdata, output_table, row.names = T, col.names = T, sep="\t", quote=F)
# Output R object
saveRDS(sc, output_rdat)

