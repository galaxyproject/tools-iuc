#!/usr/bin/env R

## Helper functions for processing variant data, especially for
## data which contains duplicate variants differing only on
## annotation.

difference.in.group <- function(lines){
    #' Find the columns containing the differences between a group of
    #' lines sharing the same POS and ALT.
    #' e.g.
    #'         CHROM  POS REF ALT  IMPACT FUNCLASS   AA   GENE
    #' 1 NC_045512.2 3037   C   T     LOW   SILENT F924 ORF1ab
    #' 2 NC_045512.2 3037   C   T     LOW   SILENT F106 ORF1ab
    #' 3 NC_045512.2 3037   C   T  MEDIUM   SILENT F106 ORF1ab
    #'
    #' should yield:
    #'
    #'            unique.name
    #' 1    3037 T (LOW|F924)
    #' 2    3037 T (LOW|F106)
    #' 3 3037 T (MEDIUM|F106)
    #'
    #' i.e. it identifies that IMPACT and AA are the differing
    #'      columns for each of the rows
    #'
    #' Ideally this function should just be used as
    #' ``tab %>% group_by(POS,ALT) %>% difference.in.group()``
    #'
    diff.colnames <- c()
    nlines = nrow(lines)
    if (nlines > 1){
        for (i in 1:(nlines-1)){
            test1 <- lines[i,]
            for (j in i:nlines){
                test2 <- lines[j,]
                diff.colnames <- c(diff.colnames,
                                   names(test1[!(test1 %in% test2)]))
            }
        }
    }
    uni.select = c("POS","ALT", diff.colnames)
    return(lines[,uni.select] %>% unite(uni.select, sep=" "))
}

splitTableAndProcess <- function(tab){
    #' Split TAB into groups sharing the same POS and ALT
    #' and create distinguishing labels.
    #'
    #' Calls the above ``difference.in.group`` for each
    #' discovered group.
    #'
    #' This function is necessary because tidyr is difficult
    #' to write custom group binding functions.
    posalts <- tab %>% group_by(POS,ALT) %>% select(POS,ALT)
    groups <- list()
    groups[[1]] <- c(1,1)
    last.pa <- paste(posalts[1,])
    for (r in 2:nrow(tab)){
        curr.pa = paste(posalts[r,])
        if (!all(last.pa == curr.pa)){
            ## end of current group, start of new
            groups[[length(groups)]][2] = r-1    ## change prev end
            groups[[length(groups)+1]] = c(r,r)  ## set start end
        }
        last.pa = curr.pa
    }
    ## split.boolean <- duplicated(tab %>%
    ##                             group_by(POS,ALT) %>%
    ##                             select(POS,ALT))
    ## split.indexes <- c(1:nrow(tab))[split.boolean]
    ## groups <- list()
    ## groups[[1]] = c(1:split.indexes[1])  ## first group
    ## for (i in 2:length(split.indexes)){
    ##     prev = split.indexes[i-1]
    ##     curr = split.indexes[i]
    ##     nindexes = c((prev+1):curr)
    ##     groups[[i]] = nindexes
    ## }
    ## num.groups =length(split.indexes)
    ## last.split = split.indexes[num.groups]
    ## last.index = nrow(tab)
    ## if (last.split < last.index){
    ##     groups[[num.groups+1]] = (last.split+1):last.index
    ## }
    as_tibble(do.call(
        "rbind",
        lapply(groups, function(g.range){
            difference.in.group(tab[unique(g.range),])
        })))
}

readAndProcess <- function(id){
    file <- (samples %>% filter(ids==id))$files
    variants <- read.table(file, header = T, sep = "\t")
    uniq.ids <- splitTableAndProcess(variants)
    stopifnot(nrow(variants) == nrow(uniq.ids))
    variants <- as_tibble(cbind(variants, uniq.ids))
    return(variants)
}
