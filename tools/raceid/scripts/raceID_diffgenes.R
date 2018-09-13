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

if (generate_extable){
                                        # differentially expressed genes in cluster

    if (!(dir.exists(outtable_dir))){
        dir.create(outtable_dir)
    }
    
    for ( n in names(cdiff) ){
        write.table(
            data.frame(
                GENEID=rownames(cdiff[[n]]),
                cdiff[[n]]
            ),
            paste(outtable_dir, "/gdiff", "_", n, ".tsv", sep=""),  # gdiff_cl.n.tsv
            row.names=FALSE,
            col.names=TRUE,
            sep="\t",
            quote=FALSE
        )        
    }
}

if (compare_clusters){   
    clust1 <- c(unlist(lapply(strsplit(clust1, "\\s*,\\s*"), as.integer)))
    clust2 <- c(unlist(lapply(strsplit(clust2, "\\s*,\\s*"), as.integer)))

    if (length(clust1) == 1){ clust1 <- clust1[[1]] }
    if (length(clust2) == 1){ clust2 <- clust2[[1]] }

    message("Performing diffgenes with cl1=%s, cl2=%s, mincount=%.0f, gene='%s'", clust1, clust2, mcount, gene_name)
    d <- diffgenes(sc,cl1 = clust1, cl2 = clust2, mincount = mcount)

    plotter("plot_diffgenes", plotdiffgenes(d,gene= gene_name))
}
