#!/usr/bin/env Rscript
# tool: pubmed_by_queries
#
# This tool uses a set of search queries to download a defined number of abstracts or
# PMIDs for search query from PubMed. PubMed's search rules and syntax apply.
#
# Input: Tab-delimited table with search queries in a column starting with "ID_",
# e.g. "ID_gene" if search queries are genes.
#
# Output: Input table with additional columns
# with PMIDs or abstracts (--abstracts) from PubMed.
#
# Usage:
# $pubmed_by_queries.R [-h] [-i INPUT] [-o OUTPUT] [-n NUMBER] [-a] [-k KEY]
#
# optional arguments:
# -h, --help                  show this help message and exit
# -i INPUT, --input INPUT     input file name. add path if file is not in working directory
# -o OUTPUT, --output OUTPUT  output file name. [default "pubmed_by_queries_output"]
# -n NUMBER, --number NUMBER  number of PMIDs or abstracts to save per ID [default "5"]
# -a, --abstract              if abstracts instead of PMIDs should be retrieved use --abstracts
# -k KEY, --key KEY           if ncbi API key is available, add it to speed up the download of PubMed data.
# For usage in Galaxy add the API key to the Galaxy user-preferences (User/ Preferences/ Manage Information).

if ("--install_packages" %in% commandArgs()) {
    print("Installing packages")
    if (!require("argparse")) install.packages("argparse", repo = "http://cran.rstudio.com/")
    if (!require("easyPubMed")) install.packages("easyPubMed", repo = "http://cran.rstudio.com/")
}

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("easyPubMed"))

parser <- ArgumentParser()
parser$add_argument("-i", "--input",
    help = "Input fie name. add path if file is not in working directory"
)
parser$add_argument("-o", "--output",
    default = "pubmed_by_queries_output",
    help = "Output file name. [default \"%(default)s\"]"
)
parser$add_argument("-n", "--number",
    type = "integer", default = 5,
    help = "Number of PMIDs (or abstracts) to save per  ID. [default \"%(default)s\"]"
)
parser$add_argument("-a", "--abstract",
    action = "store_true", default = FALSE,
    help = "If abstracts instead of PMIDs should be retrieved use --abstracts "
)
parser$add_argument("-k", "--key",
    type = "character",
    help = "If ncbi API key is available, add it to speed up the download of PubMed data. For usage in Galaxy add the API key to the Galaxy user-preferences (User/ Preferences/ Manage Information)."
)
parser$add_argument("--install_packages",
    action = "store_true", default = FALSE,
    help = "If you want to auto install missing required packages."
)
args <- parser$parse_args()

if (!is.null(args$key)) {
    if (file.exists(args$key)) {
        credentials <- read.table(args$key, quote = "\"", comment.char = "")
        args$key <- credentials[1, 1]
    }
}

max_web_tries <- 100

data <- read.delim(args$input, stringsAsFactors = FALSE)

id_col_index <- grep("ID_", names(data))


fetch_pmids <- function(data, number, pubmed_search, query, row, max_web_tries) {
    my_pubmed_url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
        "db=pubmed&retmax=", number,
        "&term=", pubmed_search$OriginalQuery,
        "&usehistory=n",
        sep = ""
    )
    # get ids
    idxml <- c()
    for (i in seq(max_web_tries)) {
        tryCatch(
            {
                id_connect <- suppressWarnings(url(my_pubmed_url, open = "rb", encoding = "UTF8"))
                idxml <- suppressWarnings(readLines(id_connect, warn = FALSE, encoding = "UTF8"))
                suppressWarnings(close(id_connect))
                break
            },
            error = function(e) {
                print(paste("Error getting URL, sleeping", 2 * i, "seconds."))
                print(e)
                Sys.sleep(time = 2 * i)
            }
        )
    }
    pmids <- c()
    for (i in seq(length(idxml))) {
        if (grepl("^<Id>", idxml[i])) {
            pmid <- custom_grep(idxml[i], tag = "Id", format = "char")
            pmids <- c(pmids, as.character(pmid[1]))
        }
    }
    if (length(pmids) > 0) {
        data[row, sapply(seq(length(pmids)), function(i) {
            paste0("PMID_", i)
        })] <- pmids
        cat(length(pmids), " PMIDs for ", query, " are added in the table.", "\n")
    }
    return(data)
}


fetch_abstracts <- function(data, number, query, pubmed_search) {
    efetch_url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
        "db=pubmed&WebEnv=", pubmed_search$WebEnv, "&query_key=", pubmed_search$QueryKey,
        "&retstart=", 0, "&retmax=", number,
        "&rettype=", "null", "&retmode=", "xml",
        sep = ""
    )
    api_key <- pubmed_search$APIkey
    if (!is.null(api_key)) {
        efetch_url <- paste(efetch_url, "&api_key=", api_key, sep = "")
    }
    # initialize
    out_data <- NULL
    try_num <- 1
    t_0 <- Sys.time()
    # Try to fetch results
    while (is.null(out_data)) {
        # Timing check: kill at 3 min
        if (try_num > 1) {
            Sys.sleep(time = 2 * try_num)
            cat(
                "Problem to receive PubMed data or error is received. Please wait. Try number:",
                try_num, "\n"
            )
        }
        t_1 <- Sys.time()
        if (as.numeric(difftime(t_1, t_0, units = "mins")) > 3) {
            message(
                "Killing the request! Something is not working. Please, try again later",
                "\n"
            )
            return(data)
        }
        # ENTREZ server connect
        out_data <- tryCatch({
            tmp_connect <- suppressWarnings(url(efetch_url,
                open = "rb",
                encoding = "UTF8"
            ))
            suppressWarnings(readLines(tmp_connect,
                warn = FALSE,
                encoding = "UTF8"
            ))
        }, error = function(e) {
            print(e)
        }, finally = {
            try(suppressWarnings(close(tmp_connect)),
                silent = TRUE
            )
        })
        # Check if error
        if (!is.null(out_data) &&
            class(out_data) == "character" &&
            grepl("<ERROR>", substr(paste(utils::head(out_data, n = 100),
                collapse = ""
            ), 1, 250))) {
            out_data <- NULL
        }
        try_num <- try_num + 1
    }
    if (is.null(out_data)) {
        message(
            "Killing the request! Something is not working. Please, try again later",
            "\n"
        )
        return(data)
    } else {
        return(out_data)
    }
}


process_xml_abstracts <- function(out_data) {
    xml_data <- paste(out_data, collapse = "")
    # articles to list
    xml_data <- strsplit(xml_data, "<PubmedArticle(>|[[:space:]]+?.*>)")[[1]][-1]
    xml_data <- sapply(xml_data, function(x) {
        # trim extra stuff at the end of the record
        if (!grepl("</PubmedArticle>$", x)) {
            x <- sub("(^.*</PubmedArticle>).*$", "\\1", x)
        }
        # Rebuid XML structure and proceed
        x <- paste("<PubmedArticle>", x)
        gsub("[[:space:]]{2,}", " ", x)
    },
    USE.NAMES = FALSE, simplify = TRUE
    )
    # titles
    titles <- sapply(xml_data, function(x) {
        x <- custom_grep(x, tag = "ArticleTitle", format = "char")
        x <- gsub("</{0,1}i>", "", x, ignore.case = T)
        x <- gsub("</{0,1}b>", "", x, ignore.case = T)
        x <- gsub("</{0,1}sub>", "", x, ignore.case = T)
        x <- gsub("</{0,1}exp>", "", x, ignore.case = T)
        if (length(x) > 1) {
            x <- paste(x, collapse = " ", sep = " ")
        } else if (length(x) < 1) {
            x <- NA
        }
        x
    },
    USE.NAMES = FALSE, simplify = TRUE
    )
    # abstracts
    abstract_text <- sapply(xml_data, function(x) {
        custom_grep(x, tag = "AbstractText", format = "char")
    },
    USE.NAMES = FALSE, simplify = TRUE
    )
    abstracts <- sapply(abstract_text, function(x) {
        if (length(x) > 1) {
            x <- paste(x, collapse = " ", sep = " ")
            x <- gsub("</{0,1}i>", "", x, ignore.case = T)
            x <- gsub("</{0,1}b>", "", x, ignore.case = T)
            x <- gsub("</{0,1}sub>", "", x, ignore.case = T)
            x <- gsub("</{0,1}exp>", "", x, ignore.case = T)
        } else if (length(x) < 1) {
            x <- NA
        } else {
            x <- gsub("</{0,1}i>", "", x, ignore.case = T)
            x <- gsub("</{0,1}b>", "", x, ignore.case = T)
            x <- gsub("</{0,1}sub>", "", x, ignore.case = T)
            x <- gsub("</{0,1}exp>", "", x, ignore.case = T)
        }
        x
    },
    USE.NAMES = FALSE, simplify = TRUE
    )
    # add title to abstracts
    if (length(titles) == length(abstracts)) {
        abstracts <- paste(titles, abstracts)
    }
    return(abstracts)
}


pubmed_data_in_table <- function(data, row, query, number, key, abstract) {
    if (is.null(query)) {
        print(data)
    }
    pubmed_search <- get_pubmed_ids(query, api_key = key)
    if (as.numeric(pubmed_search$Count) == 0) {
        cat("No PubMed result for the following query: ", query, "\n")
        return(data)
    } else if (abstract == FALSE) { # fetch PMIDs
        data <- fetch_pmids(data, number, pubmed_search, query, row, max_web_tries)
        return(data)
    } else if (abstract == TRUE) { # fetch abstracts and title text
        out_data <- fetch_abstracts(data, number, query, pubmed_search)
        abstracts <- process_xml_abstracts(out_data)
        # add abstracts to data frame
        if (length(abstracts) > 0) {
            data[row, sapply(
                seq(length(abstracts)),
                function(i) {
                    paste0("ABSTRACT_", i)
                }
            )] <- abstracts
            cat(
                length(abstracts), " abstracts for ", query, " are added in the table.",
                "\n"
            )
        }
        return(data)
    }
}

for (i in seq(nrow(data))) {
    data <- tryCatch(pubmed_data_in_table(
        data = data,
        row = i,
        query = data[i, id_col_index],
        number = args$number,
        key = args$key,
        abstract = args$abstract
    ), error = function(e) {
        print("main error")
        print(e)
        Sys.sleep(5)
    })
}

write.table(data, args$output, append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
