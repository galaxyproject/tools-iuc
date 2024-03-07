options(show.error.messages = F, error = function() {
    cat(geterrmessage(), file = stderr())
    q("no", 1, F)
})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library(qqman)
    library(optparse)
})
option_list <- list(
    make_option(c("-f", "--file"), type = "character", help = "Input file"),
    make_option("--pval",
        type = "character",
        help = "Pvalue column name", default = "P"
    ),
    make_option("--chr",
        type = "character",
        help = "Chromosome column name", default = "CHR"
    ),
    make_option("--bp",
        type = "character",
        help = "Chromosomal position column name", default = "BP"
    ),
    make_option("--snp",
        type = "character",
        help = "Snp name column name", default = "SNP"
    ),
    make_option("--name",
        type = "character",
        help = "Plot name", default = "Manhattan Plot"
    )
)
args <- parse_args(OptionParser(option_list = option_list))
file <- args$file
pvalcol <- args$pval
chrcol <- args$chr
bpcol <- args$bp
snpcol <- args$snp
name <- as.character(args$name)
data <- read.table(args$file, header = TRUE)
pdf("manhattan.pdf")
manhattan(data,
    chr = chrcol, bp = bpcol,
    p = pvalcol, snp = snpcol, main = name
)
invisible(dev.off())
