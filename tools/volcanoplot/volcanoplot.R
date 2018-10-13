#
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
  "fdr", "a", 1, "integer",
  "pvalue", "p", 1, "integer",
  "logfc", "c", 1, "integer",
  "labels", "l", 1, "integer",
  "signif", "s", 1, "double",
  "label_file", "f", 1, "character",
  "topnum", "t", 1, "integer",
  "logfc_thresh", "x", 1, "integer"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Below modified from http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html

results <- read.delim(opt$input)

results$fdr <- results[, opt$fdr]
results$Pvalue <- results[, opt$pvalue]
results$logFC <- results[, opt$logfc]
results$labels <- results[, opt$labels]

results <- mutate(results, sig=ifelse((fdr<opt$signif & logFC>opt$logfc_thresh), "Up", ifelse((fdr<opt$signif & logFC < -opt$logfc_thresh),"Down", "Not Sig")))
results <- results[order(results$Pvalue),]
if (!is.null(opt$label_file)) {
	tolabel <- read.delim(opt$label_file)
} else {
	tolabel <- head(results[results$fdr<opt$signif, ], opt$topnum)
}

pdf("out.pdf")
p <- ggplot(results, aes(logFC, -log10(Pvalue))) +
  geom_point(aes(col=sig)) +
  #geom_text_repel(data=tolabel, aes(label=labels)) +
  geom_label_repel(data=tolabel, aes(label=labels, fill=factor(sig)), colour="white", segment.colour="black", show.legend=FALSE) +
  scale_color_manual(values=c("Down"="cornflowerblue", "Not Sig"="grey", "Up"="firebrick")) +
  scale_fill_manual(values=c("Down"="cornflowerblue", "Not Sig"="grey", "Up"="firebrick")) +
  theme(panel.grid.major = element_blank(), 
  	panel.grid.minor = element_blank(),
  	panel.background = element_blank(),
  	axis.line = element_line(colour = "black"),
  	legend.key=element_blank())
print(p)
dev.off()

cat("Session information:\n\n")
sessionInfo()