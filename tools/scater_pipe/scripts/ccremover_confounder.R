#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

script_dir = args[1]
config_file = args[2]

## Load libs, common functions, source Galaxy config
source(paste(script_dir, "common.R", sep="/"))

suppressPackageStartupMessages(
    require(ccRemover)
)

## The tool only works when the mean is 0 between all cells
## so here we do just that, and add back the difference later
mean_gene_exp <- rowMeans(logcounts(sce))
t_cell_data_cen <- logcounts(sce) - mean_gene_exp

gret <- rbind(
    summary(apply(t.cell_data,1, mean)),
    summary(apply(t_cell_data_cen,1,mean))
)

rownames(gret) <- c("Before ZeroNorm:", "After ZeroNorm:")

logfile <- file("sce_log.txt", 'w')

write.table(signif(gret,6), logfile)

gene_names <- rownames(t_cell_data_cen)
cut_gene_names <- gene_names
if (regex_perform){
    cut_gene_names <- sub(regex_genenames, regex_extract, gene_names)
}

cell_cycle_gene_indices <- gene_indexer(
    cut_gene_names, species = nspecies, name_type = ntype
)

if_cc <- rep(FALSE,nrow(t_cell_data_cen)) 
if_cc[cell_cycle_gene_indices] <- TRUE

summ <- summary(if_cc)
tab <- c(as.integer(summ[[2]]), as.integer(summ[[3]]))
names(tab) <- c("False", "True")
write.table(tab, logfile, append=T)

## Begin bootstrapping process
dat <- list(x=t_cell_data_cen, if_cc=if_cc)
logs <- capture.output(
    xhat <- ccRemover(dat,
                      cutoff = ncutoff,
                      max_it = nmaxit,
                      nboot = nnboot,
                      ntop = nntop,
                      bar= F
                      )
)


## Add back the mean expression, apply the corrections to our logcounts
xhat2 <- xhat + mean_gene_exp
logcounts(sce) <- xhat2

write(logs, logfile, append=T)
saveRDS(sce, "sce_out.rds")
