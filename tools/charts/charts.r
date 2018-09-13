#!/usr/bin/Rscript

# load getopt library
library('getopt');

# convert multi parameter string (i.e. key1: value, key2: value, ...) to object
split <- function(argument){
    # process parameter string
    options <- list()
    list <- gsub("\\s","", argument)
    list <- strsplit(list, ",")
    if (length(list) > 0) {
        list <- list[[1]]
        for (entry in list) {
            pair <- strsplit(entry, ":")
            if (length(pair) > 0) {
                pair <- pair[[1]]
                if (length(pair) == 2) {
                    options[[pair[1]]] <- pair[2]
                }
            }
        }
    }
    return(options)
}

# get options, using the spec as defined by the enclosed list.
spec = matrix(c(
    'workdir',  'w', 1, 'character', 'Work directory',
    'module',   'm', 1, 'character', 'Module name',
    'input',    'i', 1, 'character', 'Input tabular file',
    'columns',  'c', 1, 'character', 'Columns string',
    'settings', 's', 1, 'character', 'Settings string',
    'output',   'o', 1, 'character', 'Output tabular file',
    'help',     'h', 0, '',          'Help',
    'verbose',  'v', 0, '',          'Verbose'
), byrow=TRUE, ncol=5);
opt = getopt(spec);

# show help
if ( !is.null(opt$help) ||
    is.null(opt$module) ||
    is.null(opt$input) ||
    is.null(opt$columns) ||
    is.null(opt$output)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1);
}

# read columns/settings
columns = split(opt$columns)
settings = split(opt$settings)

# read table
table <- read.table(opt$input, comment.char='#', fill=TRUE)

# identify module file
module_file = paste(opt$workdir, opt$module, '.r', sep='')

# source module
source(module_file)

# run module
l = wrapper (table, columns, settings)

# header
header_title <- '# title - Chart Utilities (charts)'
header_date <- paste('# date -', Sys.time(), sep=' ')
header_module <- paste('# module -', opt$module, sep=' ')
header_settings <- paste('# settings -', opt$settings, sep=' ')
header_columns <- paste('# columns -', opt$columns, sep=' ')

# check result
if (length(l) > 0) {
    # print details
    if (!is.null(opt$verbose)) {
        print ('Columns:')
        print (columns)
        print ('Settings:')
        print (settings)
        print ('Result:')
        print (l)
    }

    # create output file
    output <- file(opt$output, open='wt')
    
    # write header
    writeLines('#', output)
    writeLines(header_title, output)
    writeLines(header_date, output)
    writeLines(header_module, output)
    writeLines(header_settings, output)
    writeLines(header_columns, output)
    writeLines('#', output)
    
    # pad columns
    rows <- max(unlist(lapply(l, length)))
    padded <- lapply(l, function(col) {
        length(col) = rows;
        col
    })
    
    # write table
    write.table(padded, file=output, row.names=FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    
    # close file
    close(output)
} else {
    # print details
    print ('Columns:')
    print (columns)
    print ('Settings:')
    print (settings)
    print ('No output generated.')
}