options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library("fgsea")
  library("optparse")
})

option_list <- list(
  make_option(c("-rnk_file", "--rnk_file"), type="character", help="Path to file with differential gene expression result"),
  make_option(c("-out_tab","--out_tab"), type="character", help="Path to output file."),
  make_option(c("-gmt_file", "--gmt_file"), default="h.all.v6.1.symbols.gmt", type="character", help = "Path to Broad gmt file"),
  make_option(c("-min_size", "--min_size"), default=1, help="Minimal size of a gene set to test. All pathways below the threshold are excluded."),
  make_option(c("-max_size", "--max_size"), default=500, help="Maximal size of a gene set to test. All pathways above the threshold are excluded."),
  make_option(c("n_perm", "--n_perm"),default=1000, help="Number of permutations to do. Minimial possible nominal p-value is about 1/nperm")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
rnk_file = args$rnk_file
#category_file = args$category_file
#length_file = args$length_file
gmt_file = args$gmt_file
out_tab = args$out_tab
min_size = args$min_size
max_size = args$max_size
n_perm = args$n_perm

### Change to whatever path the gmt files are stored in
#path_to_gmt <- "."

#gmt_file <-paste(path_to_gmt, gmt_file, sep="/")

###temp values for testing

#rnk_file <- "testdata/t47d_Treatment_DEA_Prog-vs-Control_all_for_GSEA.rnk"
#gmt_file <- "h.all.v5.2.symbols.gmt"
#min_size = 5
#max_size = 500
#n_perm=1000

## Basically using the steps from the fgsea vignette

rankTab <- read.table(rnk_file,
                    header=TRUE, colClasses = c("character", "numeric"))

ranks <-rankTab[,2]
names(ranks) <- rankTab[,1]
head(ranks)

## Report an error if gmt_file not found
if(file.exists(gmt_file)) {
  pathways <- gmtPathways(gmt_file)
} else {
  cat(paste("Could not find file ", gmt_file, "in folder", path_to_gmt))
  cat("Printing contents of directory:")
  list.files()
}
fgseaRes <- fgsea(pathways, ranks, minSize=min_size, maxSize=max_size, nperm=n_perm)
write.table(fgseaRes, out_tab, sep="\t",row.names=FALSE)
