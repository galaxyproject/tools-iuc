options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
	library("argparse")
	library("gprofiler2")
	library("plyr")
})

sessionInfo()

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

set_user_agent(paste0(get_user_agent(), " galaxy"))

parser <- ArgumentParser()

parser$add_argument("--input", type="character")
parser$add_argument("--output", type="character")

parser$add_argument("--organism", type="character")
parser$add_argument("--target_ns", type="character")
parser$add_argument("--numeric_ns", type="character")

parser$add_argument("--max", type="integer")
parser$add_argument("--filter_na", action="store_true")

# Tool settings
parser$add_argument("--base_url", type="character")

args <- parser$parse_args()

query <- scan(args$input, character(), quote = "")


if (is.null(args$max)) {
	mthreshold = Inf
} else {
	mthreshold <- args$max
}

if (args$base_url != 'None') {
	set_base_url(args$base_url)
}

response <- gconvert(query
	    		  	, organism   = args$organism
					, target     = args$target_ns
					, numeric_ns = args$numeric_ns
					, mthreshold = mthreshold
					, filter_na  = args$filter_na
					)

output <- response
output.colnames = colnames(output)
write.table(output, file=args$output, quote=FALSE, sep='\t', row.names = FALSE, col.names = output.colnames)
