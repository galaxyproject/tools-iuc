# Report actual table headers, including variable sample and screening columns.
args <- commandArgs(trailingOnly = TRUE)
entries <- vapply(args, function(arg) {
    pair <- strsplit(arg, "=", fixed = TRUE)[[1]]
    header <- strsplit(readLines(pair[2], n = 1), "\t", fixed = TRUE)[[1]]
    paste0(
        '"', pair[1], '":{"metadata":{"columns":', length(header),
        ',"column_names":[', paste(encodeString(header, quote = '"'), collapse = ","), "]}}"
    )
}, character(1))
writeLines(paste0("{", paste(entries, collapse = ","), "}"), "galaxy.json")
