# wrapper
wrapper <- function(table, columns, options) {

    # initialize output list
    l <- list()

    # loop through all columns
    m <- list()
    for (key in names(columns)) {
        # load column data
        column <- as.numeric(columns[key])
        column_data <- suppressWarnings(as.numeric(as.character(table[column][[1]])))
        
        # collect vectors in list
        m <- append(m, list(column_data))
    }
    
    # identify optimal breaks
    hist_data <- hist(unlist(m), plot=FALSE)
    breaks <- hist_data$breaks;
    
    # add as first column
    l <- append(l, list(breaks[2: length(breaks)]))
    
    # loop through all columns
    for (key in seq(m)) {
        # load column data
        column_data <- m[[key]]
        
        # create hist data
        hist_data <- hist(column_data, breaks=breaks, plot=FALSE)
        
        # normalize densities
        count_sum <- sum(hist_data$counts)
        if (count_sum > 0) {
            hist_data$counts = hist_data$counts / count_sum
        }
        
        # collect vectors in list
        l <- append(l, list(hist_data$counts))
    }
    
    # return
    return (l)
}
