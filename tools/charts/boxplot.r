wrapper <- function(table, columns, options) {

    # initialize output list
    l <- list()

    # loop through all columns
    for (key in names(columns)) {
        # load column data
        column <- as.numeric(columns[key])
        column_data <- suppressWarnings(as.numeric(as.character(table[column][[1]])))
    
        # create hist data
        data <- boxplot(column_data, plot=FALSE)
        
        # collect vectors in list
        l <- append(l, list(data$stats))
    }
    
    # return
    return (l)
}
