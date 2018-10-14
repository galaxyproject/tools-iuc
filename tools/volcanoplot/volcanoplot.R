# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(dplyr)
    library(getopt)
    library(ggplot2)
    library(ggrepel)
})

options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

spec <- matrix(c(
    "input", "i", 1, "character",
    "fdr_col", "a", 1, "integer",
    "pval_col", "p", 1, "integer",
    "lfc_col", "c", 1, "integer",
    "label_col", "l", 1, "integer",
    "signif_thresh", "s", 1, "double",
    "lfc_thresh", "x", 1, "double",
    "label_file", "f", 1, "character",
    "top_num", "t", 1, "integer",
    "title", "T", 1, "character",
    "xlab", "X", 1, "character",
    "ylab", "Y", 1, "character",
    "legend", "L", 1, "character"),
    byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Below modified from http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html

results <- read.delim(opt$input)
results$fdr <- results[, opt$fdr_col]
results$Pvalue <- results[, opt$pval_col]
results$logFC <- results[, opt$lfc_col]
results$labels <- results[, opt$label_col]

results <- mutate(results, sig=ifelse((fdr<opt$signif_thresh & logFC>opt$lfc_thresh), "Up", ifelse((fdr<opt$signif_thresh & logFC < -opt$lfc_thresh),"Down", "Not Sig")))
results <- results[order(results$Pvalue),]
if (!is.null(opt$label_file)) {
    labelfile <- read.delim(opt$label_file)
    tolabel <- filter(results, labels %in% labelfile[, 1])
} else if (!is.null(opt$top_num)) {
    tolabel <- filter(results, fdr<opt$signif_thresh) %>% top_n(opt$top_num)
} else {
    tolabel <- filter(results, fdr<opt$signif_thresh)
}

pdf("out.pdf")
p <- ggplot(results, aes(logFC, -log10(Pvalue))) +
    geom_point(aes(col=sig)) +
    geom_label_repel(data=tolabel, aes(label=labels, fill=factor(sig)), colour="white", segment.colour="black", show.legend=FALSE) +
    scale_color_manual(values=c("Down"="cornflowerblue", "Not Sig"="grey", "Up"="firebrick")) +
    scale_fill_manual(values=c("Down"="cornflowerblue", "Not Sig"="grey", "Up"="firebrick")) +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key=element_blank())
    if (!is.null(opt$title)) {
        p <- p + ggtitle(opt$title)
    }
    if (!is.null(opt$xlab)) {
        p <- p + xlab(opt$xlab)
    }
    if (!is.null(opt$ylab)) {
        p <- p + ylab(opt$ylab)
    }
    if (!is.null(opt$legend)) {
        p <- p + labs(colour=opt$legend)
    } else {
        p <- p + labs(colour="")
    }
print(p)
dev.off()

cat("Session information:\n\n")
sessionInfo()