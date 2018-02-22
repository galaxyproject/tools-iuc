#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

                                        # Load libs, common functions, source RaceID,
                                        # Galaxy Params, and read input data (sc)
source(paste(script_dir, "common.R", sep="/"))

                                        # Read input data
message("Count matrix with %.0f cells and %.0f genes", dim(sc@fdata)[1], dim(sc@fdata)[2])

cdiff <- clustdiffgenes(sc, pvalue=c_pval)

if (generate_finheat){
    png("plot_finalheat.png")
    x <- clustheatmap(sc,final=TRUE,hmethod="single")
    dev.off()
}

if (generate_extable){
                                        # differentially expressed genes in cluster
    for ( n in names(cdiff) ){
        write.table(data.frame(
            GENEID=rownames(cdiff[[n]]),cdiff[[n]]),
            paste(paste("gdiff",
                        sub("\.","\\_",n),
                        sep="_"),".xls",
                  sep=""),
            row.names=FALSE,
            col.names=TRUE,
            sep="\t",quote=FALSE
            )
    }
}

if (compare_clusters){
    clust1 <- c(unlist(lapply(strsplit(clust1, "\\s*,\\s*"), as.integer)))
    clust2 <- c(unlist(lapply(strsplit(clust2, "\\s*,\\s*"), as.integer)))

    if (length(clust1) == 1){ clust1 <- clust1[[1]] }
    if (length(clust2) == 1){ clust2 <- clust2[[1]] }
    
    d <- diffgenes(sc,cl1 = clust1, cl2 = clust2, mincount = mcount)

    png("plot_diffgenes.png")
    plotdiffgenes(d,gene=gene_name)
    dev.off()
}


if (generate_robject){
    message("Saving SC object")
    saveRDS(sc, output_rdat)
}
