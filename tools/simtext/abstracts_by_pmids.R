#!/usr/bin/env Rscript
# TOOL2 abstracts_by_pmids
#
# This tool retrieves for all PMIDs in each row of a table the according abstracts and saves them in additional columns.
#
# Input: Tab-delimited table with columns containing PMIDs. The names of the PMID columns should start with “PMID”, e.g. “PMID_1”, “PMID_2” etc.
#
# Output: Input table with additional columns containing abstracts corresponding to the PMIDs from PubMed.
# The abstract columns are called "ABSTRACT_1", "ABSTARCT_2" etc.
#
# Usage: $ T2_abstracts_by_pmid.R [-h] [-i INPUT] [-o OUTPUT]
#
# optional arguments:
# -h, --help                 show help message
# -i INPUT, --input INPUT    input file name. add path if file is not in working directory
# -o OUTPUT, --output OUTPUT output file name. [default "T2_output"]


if ("--install_packages" %in% commandArgs()) {
    print("Installing packages")
    if (!require("argparse")) install.packages("argparse", repo = "http://cran.rstudio.com/")
    if (!require("reutils")) install.packages("reutils", repo = "http://cran.rstudio.com/")
    if (!require("easyPubMed")) install.packages("easyPubMed", repo = "http://cran.rstudio.com/")
    if (!require("textclean")) install.packages("textclean", repo = "http://cran.rstudio.com/")
}

suppressPackageStartupMessages(library("argparse"))
library("reutils")
suppressPackageStartupMessages(library("easyPubMed"))
suppressPackageStartupMessages(library("textclean"))

parser <- ArgumentParser()
parser$add_argument("-i", "--input",
    help = "input fie name. add path if file is not in workind directory"
)
parser$add_argument("-o", "--output",
    default = "abstracts_by_pmids_output",
    help = "output file name. [default \"%(default)s\"]"
)
parser$add_argument("--install_packages",
    action = "store_true", default = FALSE,
    help = "If you want to auto install missing required packages."
)

args <- parser$parse_args()

data <- read.delim(args$input, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
pmids_cols_index <- grep("PMID", names(data))

fetch_abstracts <- function(pmids, row) {
    efetch_result <- NULL
    try_num <- 1
    t_0 <- Sys.time()

    while (is.null(efetch_result)) {
        # Timing check: kill at 3 min
        if (try_num > 1) {
            Sys.sleep(time = 1 * try_num)
            cat("Problem to receive PubMed data or error is received. Please wait. Try number: ", try_num, "\n")
        }

        t_1 <- Sys.time()

        if (as.numeric(difftime(t_1, t_0, units = "mins")) > 3) {
            message("Killing the request! Something is not working. Please, try again later", "\n")
            return(data)
        }

        efetch_result <- tryCatch(
            {
                suppressWarnings(efetch(uid = pmids, db = "pubmed", retmode = "xml"))
            },
            error = function(e) {
                NULL
            }
        )

        if (!is.null(as.list(efetch_result$errors)$error)) {
            if (as.list(efetch_result$errors)$error == "HTTP error: Status 400; Bad Request") {
                efetch_result <- NULL
            }
        }

        try_num <- try_num + 1
    } # while loop end

    # articles to list
    xml_data <- strsplit(efetch_result$content, "<PubmedArticle(>|[[:space:]]+?.*>)")[[1]][-1]
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

    abstracts <- as.character(abstracts)

    if (length(abstracts) > 0) {
        data[row, sapply(seq(length(abstracts)), function(i) {
            paste0("ABSTRACT_", i)
        })] <- abstracts
        cat(length(abstracts), " abstracts for PMIDs of row ", row, " are added in the table.", "\n")
    }

    return(data)
}


for (row in seq(nrow(data))) {
    pmids <- as.character(unique(data[row, pmids_cols_index]))
    pmids <- pmids[!pmids == "NA"]

    if (length(pmids) > 0) {
        data <- tryCatch(fetch_abstracts(pmids, row),
            error = function(e) {
                Sys.sleep(3)
            }
        )
    } else {
        print(paste("No PMIDs in row", row))
    }
}
write.table(data, args$output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
