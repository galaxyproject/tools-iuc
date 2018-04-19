# load sparse matrix package
suppressPackageStartupMessages(library('Matrix'))

# access a numeric column
get_numeric <- function(table, column_key) {
    column <- as.numeric(column_key)
    column_data <- suppressWarnings(as.numeric(as.character(table[column][[1]])))
    return (c(column_data))
}

# access a label column
get_label <- function(table, column_key) {
    column <- as.numeric(column_key)
    column_data <- as.character(table[column][[1]])
    return (c(column_data))
}

# inflate three columns into matrix
matrify <- function (data) {
    if (ncol(data) != 3)
        stop('Data frame must have three column format')
    plt <- data[, 1]
    spc <- data[, 2]
    abu <- data[, 3]
    plt.codes <- levels(factor(plt))
    spc.codes <- levels(factor(spc))
    taxa <- Matrix(0, nrow=length(plt.codes), ncol=length(spc.codes), sparse=TRUE)
    row <- match(plt, plt.codes)
    col <- match(spc, spc.codes)
    for (i in 1:length(abu)) {
        taxa[row[i], col[i]] <- abu[i]
    }
    colnames(taxa) <- spc.codes
    rownames(taxa) <- plt.codes
    taxa
}

# flatten data.frame into three column format
flatten <- function(my_matrix) {
    summ <-summary(my_matrix)
    summ <- data.frame(i=rownames(my_matrix)[summ$i], j=colnames(my_matrix)[summ$j], x=summ$x)
    summ
}

# wrapper
wrapper <- function(table, columns, options) {

    # initialize output list
    l <- list()

    # get number of columns
    n = length(columns)
    
    # consistency check
    if (n %% 3 != 0) {
        print ('heatmap::wrapper() - Data not consistent (n mod 3 != 0)')
        return (l)
    }
    
    # create index sequence
    index = seq(1, n, by=3)
    
    # get keys
    keys = names(columns)
    
    # loop through blocks
    for (i in index) {
        # create columns
        ci <- get_label(table, columns[keys[i]])
        cj <- get_label(table, columns[keys[i+1]])
        cx <- get_numeric(table, columns[keys[i+2]])
        
        # create a frame from columns
        my_frame <- data.frame(ci=ci, cj=cj, cx=cx)
        
        # create matrix out of the frame
        my_matrix <- matrify(my_frame)
        
        # create/cluster matrix
        row_order <- hclust(dist(my_matrix))$order
        col_order <- hclust(dist(t(my_matrix)))$order
        
        # reorder matrix
        my_matrix <- my_matrix[row_order, col_order]
        
        # transform back to three columns
        my_flatmatrix = flatten(my_matrix)
        
        # append to result list
        l <- append(l, list(my_flatmatrix$i))
        l <- append(l, list(my_flatmatrix$j))
        l <- append(l, list(my_flatmatrix$x))
    }
    
    # return
    return (l)
}
