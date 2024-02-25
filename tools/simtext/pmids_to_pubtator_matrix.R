#!/usr/bin/env Rscript
# tool: pmids_to_pubtator_matrix
#
# The tool uses all PMIDs per row and extracts "Gene", "Disease", "Mutation", "Chemical" and "Species" terms of the
# corresponding abstracts, using PubTator annotations. The user can choose from which categories terms should be extracted.
# The extracted terms are united in one large binary matrix, with 0= term not present in abstracts of that row and 1= term
# present in abstracts of that row. The user can decide if the extracted scientific terms should be extracted and used as
# they are or if they should be grouped by their geneIDs/ meshIDs (several terms can often be grouped into one ID).
# äAlso, by default all terms are extracted, otherwise the user can specify a number of most frequent words to be extracted per row.
#
# Input: Output of abstracts_by_pmids or tab-delimited table with columns containing PMIDs.
# The names of the PMID columns should start with "PMID", e.g. "PMID_1", "PMID_2" etc.
#
# Output: Binary matrix in that each column represents one of the extracted terms.
#
# usage: $ pmids_to_pubtator_matrix.R [-h] [-i INPUT] [-o OUTPUT] [-n NUMBER]
# [-c {Genes,Diseases,Mutations,Chemicals,Species} [{Genes,Diseases,Mutations,Chemicals,Species} ...]]
#
# optional arguments:
#   -h, --help                 show help message
#   -i INPUT, --input INPUT    input file name. add path if file is not in workind directory
#   -n NUMBER, --number NUMBER Number of most frequent terms/IDs to extract. By default all terms/IDs are extracted.
#   -o OUTPUT, --output OUTPUT output file name. [default "pmids_to_pubtator_matrix_output"]
#   -c {Gene,Disease,Mutation,Chemical,Species} [{Genes,Diseases,Mutations,Chemicals,Species} ...], --categories {Gene,Disease,Mutation,Chemical,Species} [{Gene,Disease,Mutation,Chemical,Species} ...]
#      Pubtator categories that should be considered.  [default "('Gene', 'Disease', 'Mutation','Chemical')"]

if ("--install_packages" %in% commandArgs()) {
    print("Installing packages")
    if (!require("argparse")) install.packages("argparse", repo = "http://cran.rstudio.com/")
    if (!require("stringr")) install.packages("stringr", repo = "http://cran.rstudio.com/")
    if (!require("RCurl")) install.packages("RCurl", repo = "http://cran.rstudio.com/")
    if (!require("stringi")) install.packages("stringi", repo = "http://cran.rstudio.com/")
}

suppressPackageStartupMessages(library("argparse"))
library("stringr")
library("RCurl")
library("stringi")

parser <- ArgumentParser()

parser$add_argument("-i", "--input",
    help = "input fie name. add path if file is not in workind directory"
)
parser$add_argument("-o", "--output",
    default = "pmids_to_pubtator_matrix_output",
    help = "output file name. [default \"%(default)s\"]"
)
parser$add_argument("-c", "--categories",
    choices = c("Gene", "Disease", "Mutation", "Chemical", "Species"), nargs = "+",
    default = c("Gene", "Disease", "Mutation", "Chemical"),
    help = "Pubtator categories that should be considered. [default \"%(default)s\"]"
)
parser$add_argument("-b", "--byid",
    action = "store_true", default = FALSE,
    help = "If you want to find common gene IDs / mesh IDs instead of scientific terms."
)
parser$add_argument("-n", "--number",
    default = NULL, type = "integer",
    help = "Number of most frequent terms/IDs to extract. By default all terms/IDs are extracted."
)
parser$add_argument("--install_packages",
    action = "store_true", default = FALSE,
    help = "If you want to auto install missing required packages."
)

args <- parser$parse_args()


data <- read.delim(args$input, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

pmid_cols_index <- grep(c("PMID"), names(data))
word_matrix <- data.frame()
dict_table <- data.frame()
pmids_count <- 0
pubtator_max_ids <- 100


merge_pubtator_table <- function(out_data, table) {
    out_data <- unlist(strsplit(out_data, "\n", fixed = T))
    for (i in 3:length(out_data)) {
        temps <- unlist(strsplit(out_data[i], "\t", fixed = T))
        if (length(temps) == 5) {
            temps <- c(temps, NA)
        }
        if (length(temps) == 6) {
            table <- rbind(table, temps)
        }
    }
    return(table)
}


get_pubtator_terms <- function(pmids) {
    table <- NULL
    for (pmid_split in split(pmids, ceiling(seq_along(pmids) / pubtator_max_ids))) {
        out_data <- NULL
        try_num <- 1
        t_0 <- Sys.time()
        while (TRUE) {
            # Timing check: kill at 3 min
            if (try_num > 1) {
                cat("Connection problem. Please wait. Try number:", try_num, "\n")
                Sys.sleep(time = 2 * try_num)
            }
            try_num <- try_num + 1
            t_1 <- Sys.time()
            if (as.numeric(difftime(t_1, t_0, units = "mins")) > 3) {
                message("Killing the request! Something is not working. Please, try again later", "\n")
                return(table)
            }
            out_data <- tryCatch({
                getURL(paste("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=",
                    paste(pmid_split, collapse = ","),
                    sep = ""
                ))
            }, error = function(e) {
                print(e)
                next
            }, finally = {
                Sys.sleep(0)
            })
            if (!is.null(out_data)) {
                table <- merge_pubtator_table(out_data, table)
                break
            }
        }
    }
    return(table)
}

extract_category_terms <- function(table, categories) {
    index_categories <- c()
    categories <- as.character(unlist(categories))
    if (ncol(table) == 6) {
        for (i in categories) {
            tmp_index <- grep(TRUE, i == as.character(table[, 5]))
            if (length(tmp_index) > 0) {
                index_categories <- c(index_categories, tmp_index)
            }
        }
        table <- as.data.frame(table, stringsAsFactors = FALSE)
        table <- table[index_categories, c(4, 6)]
        table <- table[!is.na(table[, 2]), ]
        table <- table[!(table[, 2] == "NA"), ]
        table <- table[!(table[, 1] == "NA"), ]
    } else {
        return(NULL)
    }
}

extract_frequent_ids_or_terms <- function(table) {
    if (is.null(table)) {
        return(NULL)
        break
    }
    if (args$byid) {
        if (!is.null(args$number)) {
            # retrieve top X mesh_ids
            table_mesh <- as.data.frame(table(table[, 2]))
            colnames(table_mesh)[1] <- "mesh_id"
            table <- table[order(table_mesh$Freq, decreasing = TRUE), ]
            table <- table[1:min(args$number, nrow(table_mesh)), ]
            table_mesh$mesh_id <- as.character(table_mesh$mesh_id)
            # subset table for top X mesh_ids
            table <- table[which(as.character(table$V6) %in% as.character(table_mesh$mesh_id)), ]
            table <- table[!duplicated(table[, 2]), ]
        } else {
            table <- table[!duplicated(table[, 2]), ]
        }
    } else {
        if (!is.null(args$number)) {
            table[, 1] <- tolower(as.character(table[, 1]))
            table <- as.data.frame(table(table[, 1]))
            colnames(table)[1] <- "term"
            table <- table[order(table$Freq, decreasing = TRUE), ]
            table <- table[1:min(args$number, nrow(table)), ]
            table$term <- as.character(table$term)
        } else {
            table[, 1] <- tolower(as.character(table[, 1]))
            table <- table[!duplicated(table[, 1]), ]
        }
    }
    return(table)
}


# for all PMIDs of a row get PubTator terms and add them to the matrix
for (i in seq(nrow(data))) {
    pmids <- as.character(data[i, pmid_cols_index])
    pmids <- pmids[!pmids == "NA"]
    if (pmids_count > 10000) {
        cat("Break (10s) to avoid killing of requests. Please wait.", "\n")
        Sys.sleep(10)
        pmids_count <- 0
    }
    pmids_count <- pmids_count + length(pmids)
    # get puptator terms and process them with functions
    if (length(pmids) > 0) {
        table <- get_pubtator_terms(pmids)
        table <- extract_category_terms(table, args$categories)
        table <- extract_frequent_ids_or_terms(table)
        if (!is.null(table)) {
            colnames(table) <- c("term", "mesh_id")
            # add data in binary matrix
            if (args$byid) {
                mesh_ids <- as.character(table$mesh_id)
                if (length(mesh_ids) > 0) {
                    word_matrix[i, mesh_ids] <- 1
                    cat(length(mesh_ids), " IDs for PMIDs of row", i, " were added", "\n")
                    # add data in dictionary
                    dict_table <- rbind(dict_table, table)
                    dict_table <- dict_table[!duplicated(as.character(dict_table[, 2])), ]
                }
            } else {
                terms <- as.character(table[, 1])
                if (length(terms) > 0) {
                    word_matrix[i, terms] <- 1
                    cat(length(terms), " terms for PMIDs of row", i, " were added.", "\n")
                }
            }
        }
    } else {
        cat("No terms for PMIDs of row", i, " were found.", "\n")
    }
}

if (args$byid) {
    # change column names of matrix: exchange mesh ids/ids with term
    index_names <- match(names(word_matrix), as.character(dict_table[[2]]))
    names(word_matrix) <- dict_table[index_names, 1]
}

colnames(word_matrix) <- gsub("[^[:print:]]", "", colnames(word_matrix))
colnames(word_matrix) <- gsub('\"', "", colnames(word_matrix), fixed = TRUE)

# merge duplicated columns
word_matrix <- as.data.frame(do.call(cbind, by(t(word_matrix), INDICES = names(word_matrix), FUN = colSums)))

# save binary matrix
word_matrix <- as.matrix(word_matrix)
word_matrix[is.na(word_matrix)] <- 0
cat("Matrix with ", nrow(word_matrix), " rows and ", ncol(word_matrix), " columns generated.", "\n")
write.table(word_matrix, args$output, row.names = FALSE, sep = "\t", quote = FALSE)
