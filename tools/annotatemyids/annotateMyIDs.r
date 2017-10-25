options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-id_file", "--id_file"), type="character", help="Path to file with IDs to convert"),
  make_option(c("-out_tab","--out_tab"), type="character", help="Path to output file."),
  make_option(c("-id_type","--id_type"),type="character",help="Type of the incoming IDs"),
  make_option(c("-organism","--organism"), type="character",help="Which organism the IDs belong to"),
  make_option(c("-include_go","--include_go"),type="logical",default=TRUE,help="if TRUE, include GO IDs in the output"),
  make_option(c("-include_kegg","--include_kegg"),type="logical",default=TRUE,help="If TRUE, include KEGG pathways in the output"),
  make_option(c("-file_has_header","--file_has_header"),type="logical",default=TRUE,help="If this option is set to TRUE, the tool will assume that the ranked gene-list has a column heading and the gene names commence on the second line")
  
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
id_file <- args$id_file
out_tab <- args$out_tab
id_type <- args$id_type
organism <- args$organism
include_go <- args$include_go
include_kegg <- args$include_kegg
file_has_header = args$file_has_header

ids <- as.character(read.table(id_file)[,1],header=file_has_header)

if(organism == "hs"){
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  db <- org.Hs.eg.db
} else if (organism == "mm"){
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  db <- org.Mm.eg.db
} else if (organism == "dm"){
  suppressPackageStartupMessages(library(org.Dm.eg.db))
  db <- org.Dm.eg.db
} else if (organism == "dr"){
  suppressPackageStartupMessages(library(org.Dr.eg.db))
  db <- org.Dr.eg.db
  } else {
    cat(paste("Organism type not supported", organism))
  }

columns <- c("ENSEMBL","ENTREZID","SYMBOL","GENENAME")
if (include_go) columns <- c(columns, "GO")
if (include_kegg) columns <- c(columns, "PATH")

result <- select(db, keys = ids,keytype = id_type,columns = columns )

write.table(result, file=out_tab,sep="\t",row.names=FALSE,quote=FALSE)
