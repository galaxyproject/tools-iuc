# wrapper
wrapper <- function(table, columns, options) {

    # initialize output list
    l <- list()

    # loop through all columns
    m <- list()
    for (key in names(columns)) {
        # load column data
        column <- as.numeric(columns[key])
        
        # ensure string column
        column_data <- as.character(table[column][[1]])
    
        # collect vectors in list
        m <- append(m, list(column_data))
    }
    
    # get alphabetically sorted bins
    bins <- sort(unique(unlist(m)))
    
    # add first column
    l <- append(l, list(bins))
    
    # loop through all columns
    for (key in seq(m)) {
        # reset bins
        bins = sapply(bins, function(v) { 0 })
        
        # load column data
        column_data <- m[[key]]
        
        # create hist data
        table_data <- table(column_data)
        
        # transfer counts to bins
        for (id in names(table_data)) {
            bins[id] <- table_data[id]
        }
        
        # normalize densities
        total <- length(column_data)
        if (total > 0) {
            bins = bins / total
        }
        
        # collect vectors in list
        l <- append(l, list(bins))
    }

    # return
    return (l)
}
