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
    "legend", "L", 1, "character",
    "llabs", "z", 1, "character"),
    byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Below modified from http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html

results <- read.delim(opt$input)
results$fdr <- results[, opt$fdr_col]
results$Pvalue <- results[, opt$pval_col]
results$logFC <- results[, opt$lfc_col]
results$labels <- results[, opt$label_col]
label_down <- unlist(strsplit(opt$llabs, split=","))[1]
label_notsig <- unlist(strsplit(opt$llabs, split=","))[2]
label_up <- unlist(strsplit(opt$llabs, split=","))[3]
colours <- setNames(c("cornflowerblue","grey","firebrick"),c(label_down,label_notsig,label_up))

results <- mutate(results, sig=ifelse((fdr<opt$signif_thresh & logFC>opt$lfc_thresh), label_up, ifelse((fdr<opt$signif_thresh & logFC < -opt$lfc_thresh),label_down, label_notsig)))
results <- results[order(results$Pvalue),]
if (!is.null(opt$label_file)) {
    labelfile <- read.delim(opt$label_file)
    # label genes specified in file
    tolabel <- filter(results, labels %in% labelfile[, 1])
} else if (is.null(opt$top_num)) {
    # label all significant genes
    tolabel <- filter(results, sig != label_notsig)
} else if (opt$top_num > 0) {
    # label only top significant genes
    tolabel <- filter(results, sig != label_notsig) %>%
    top_n(n=-opt$top_num, Pvalue)
} else if (opt$top_num == 0) {
    # no labels
    tolabel <- NULL
}

pdf("out.pdf")
p <- ggplot(results, aes(logFC, -log10(Pvalue))) +
    geom_point(aes(col=sig)) +
    scale_color_manual(values=colours) +
    scale_fill_manual(values=colours) +
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
if (!is.null(tolabel)) {
    p <- p + geom_label_repel(data=tolabel, aes(label=labels, fill=factor(sig)), colour="white", segment.colour="black", show.legend=FALSE)
}

print(p)
dev.off()

cat("Session information:\n\n")
sessionInfo()