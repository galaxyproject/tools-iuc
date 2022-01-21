options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
    library("pathview")
    library("optparse")
})

sessionInfo()

option_list <- list(
  make_option(c("--pathway_id"), type="character", default=NULL, help="Path to tabular file with gene data"),
  make_option(c("--pathway_id_fp"), type="character", default=NULL, help="Path to tabular file with pathway ids"),
  make_option(c("--pathway_id_header"), type="logical", default=FALSE, help="Header for tabular file with pathway ids"),
  make_option(c("--species"), type="character", default="hsa", help="KEGG code, scientific name or the common name of the species"),
  make_option(c("--gene_data"), type="character", default=NULL, help="Path to tabular file with gene data"),
  make_option(c("--gd_header"), type="logical", default=FALSE, help="Header for the gene data file"),
  make_option(c("--gene_idtype"), type="character", default="entrez", help="ID type used for the gene data"),
  make_option(c("--cpd_data"), type="character", default=NULL, help="Path to tabular file with compound data"),
  make_option(c("--cpd_header"), type="logical", default=FALSE, help="Header for the compound data file"),
  make_option(c("--cpd_idtype"), type="character", default="kegg", help="ID type used for the compound data"),
  make_option(c("--multi_state"), type="logical", default=TRUE, help="Are the gene and compound data paired?"),
  make_option(c("--match_data"), type="logical", default=TRUE, help="Are the gene and compound data paired?"),
  make_option(c("--kegg_native"), type="logical", default=TRUE, help="Render pathway graph as native KEGG grap? Alternative is the Graphviz layout"),
  make_option(c("--same_layer"), type="logical", default=TRUE, help="Plot on same layer?"),
  make_option(c("--map_null"), type="logical", default=TRUE, help="Map the NULL gene or compound data to pathway?"),
  make_option(c("--split_group"), type="logical", default=FALSE, help="Split node groups into individual nodes?"),
  make_option(c("--expand_node"), type="logical", default=FALSE, help="Expand multiple-gene nodes into single-gene nodes?"),
  make_option(c("--sign_pos"), type="character", default="bottomright", help="Position of pathview signature")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)
print(args)

read_table = function(fp, header, rownames=1, colclasses=NA){
  table = read.table(fp, header=header, sep="\t", row.names=rownames, colClasses=colclasses)
  # transform to vector if only one column
  if(dim(table)[2] == 1){
    names = rownames(table)
    table = table[,1]
    names(table) = names
  }
  return(table)
}

get_table = function(fp, header){
  table = NULL
  if(!is.null(fp)){
    table = read_table(fp, header, rownames=1)
  }
  return(table)
}

# load gene_data file
gene_data = get_table(args$gene_data, args$gd_header)

# load compound data file
cpd_data = get_table(args$cpd_data, args$cpd_header)

run_pathview = function(pathway_id){
  pathview(
    pathway.id=pathway_id,
    gene.data=gene_data,
    gene.idtype=args$gene_idtype,
    cpd.data=cpd_data,
    cpd.idtype=args$cpd_idtype,
    species=args$species,
    multi.state=args$multi_state,
    match.data=args$match_data,
    kegg.native=args$kegg_native,
    same.layer=args$same_layer,
    split.group=args$split_group,
    expand.node=args$expand_node,
    sign.pos=args$sign_pos,
    map.null=args$map_null)
}

# get pathway ids
if(!is.null(args$pathway_id)){
  run_pathview(args$pathway_id)
} else {
  pthws = read_table(args$pathway_id_fp, args$pathway_id_header, rownames=NULL, colclasses="character")
  for(p in pthws){
    run_pathview(p)
  }
}
