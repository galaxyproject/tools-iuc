options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
	library("argparse")
	library("gprofiler2")
	library("ggplot2")
})

sessionInfo()

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

set_user_agent(paste0(get_user_agent(), " galaxy"))

parser <- ArgumentParser()

parser$add_argument("--input", type="character")
parser$add_argument("--output", type="character")

parser$add_argument("--organism", type="character")
parser$add_argument("--ordered_query", action="store_true")

parser$add_argument("-p", "--plot", type="character")

# Advanced options
parser$add_argument("--significant", action="store_true")
parser$add_argument("--exclude_iea", action="store_true")
parser$add_argument("--measure_underrepresentation", action="store_true")
parser$add_argument("--evcodes", action="store_true")
parser$add_argument("--user_threshold", type="double")
parser$add_argument("--correction_method", type="character")
parser$add_argument("--domain_scope", type="character")
parser$add_argument("--custom_bg", type="character")
parser$add_argument("--numeric_ns", type="character")

# Datasources
parser$add_argument("--sources", type="character", action="append")

# Tool settings
parser$add_argument("--base_url", type="character")

args <- parser$parse_args()

query <- scan(args$input, character(), quote = "")

if (args$custom_bg != 'None') {
	custom_bg <- scan(args$custom_bg, character(), quote = "")
} else {
	custom_bg = NULL
}

if (args$base_url != 'None') {
	set_base_url(args$base_url)
}


if (length(args$sources) > 0) {
	sources <- unlist(args$sources)
}

response <- gost(query
				, organism                    = args$organism
				, ordered_query               = args$ordered_query
				, significant                 = args$significant
				, exclude_iea                 = args$exclude_iea
				, measure_underrepresentation = args$measure_underrepresentation
				, evcodes                     = args$evcodes
				, user_threshold              = args$user_threshold
				, correction_method           = args$correction_method
				, domain_scope                = args$domain_scope
				, custom_bg                   = custom_bg
				, numeric_ns                  = args$numeric_ns
				, sources                     = sources
				)

output <- response$result
output$parents <- vapply(output$parents, paste, collapse = ",", character(1L))
output.colnames = colnames(output)
write.table(output, file=args$output, quote=FALSE, sep='\t', row.names = FALSE, col.names = output.colnames)

if (args$plot != 'None') {
	image.name <- strsplit(args$plot, "\\.")
	plot <- gostplot(response, interactive = FALSE)
	ggsave(plot, filename = paste0(unlist(image.name)[1], ".png"))

	# publish_gostplot(plot, highlight_terms = NULL, filename = paste0(unlist(image.name)[1], ".png"))
	file.rename(paste0(unlist(image.name)[1], ".png"), paste0(unlist(image.name)[1], ".dat"))
}
