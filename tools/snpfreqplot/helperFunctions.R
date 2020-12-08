#!/usr/bin/env R

## Helper functions for processing variant data, especially for
## data which contains duplicate variants differing only on
## annotation.

difference_in_group <- function(lines) {
    #' Find the columns containing the differences between a
    #' group of lines sharing the same POS and ALT.
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
    #' ``tab %>% group_by(POS, ALT) %>% difference_in_group()``
    #'
    diff.colnames <- c()
    nlines <- nrow(lines)
    if (nlines > 1) {
        for (i in 1:(nlines - 1)) {
            test1 <- lines[i, ]
            for (j in i:nlines) {
                test2 <- lines[j, ]
                diff.colnames <- c(diff.colnames,
                                   names(test1[!(test1 %in% test2)]))
            }
        }
    }
    group_select <- c("POS", "REF", "ALT", diff.colnames)
    return(lines[, group_select] %>% unite(group_select, sep = " ")) # nolint
}

split_table_and_process <- function(tab) {
    #' Split TAB into groups sharing the same POS and ALT
    #' and create distinguishing labels.
    #'
    #' Calls the above ``difference_in_group`` for each
    #' discovered group.
    #'
    #' This function is necessary because tidyr is difficult
    #' to write custom group binding functions.
    group_ind <- tab %>% group_by(POS, REF, ALT) %>% select(POS, REF, ALT) # nolint
    nlines <- nrow(tab)
    groups <- list()
    groups[[1]] <- c(1, 1)
    last_pa <- paste(group_ind[1, ])
    for (r in 2:nlines) {
        curr_pa <- paste(group_ind[r, ])
        group_ind_diff_between_lines <- !all(last_pa == curr_pa)
        if (group_ind_diff_between_lines) {
            ## end of current group, start of new
            groups[[length(groups)]][2] <- r - 1     ## change prev end
            groups[[length(groups) + 1]] <- c(r, r)  ## set (start, end)
        } else if (r == nlines) {
            ## i.e. if the very last line shares
            ## the same POS REF ALT as the one before,
            ## close current group.
            groups[[length(groups)]][2] <- r
        }
        last_pa <- curr_pa
    }
    as_tibble(do.call(
        "rbind",
        lapply(groups, function(grange) {
            expand_range <- grange[1]:grange[2]
            difference_in_group(tab[expand_range, ])
        })
    ))
}

read_and_process <- function(id) {
    file <- (samples %>% filter(ids == id))$files    # nolint
    variants <- read.table(file, header = T, sep = "\t")
    uniq_ids <- split_table_and_process(variants)
    if (nrow(variants) != nrow(uniq_ids)) {
        stop(paste0(id, " '", file, "' failed: ", file, "\"",
                    "nrow(variants)=", nrow(variants),
                    " but nrow(uniq_ids)=", nrow(uniq_ids)))
    }
    variants <- as_tibble(cbind(variants, uniq_ids)) # nolint
    return(variants)
}
