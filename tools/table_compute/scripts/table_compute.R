#!/usr/bin/env R
VERSION = "0.5"

args = commandArgs(trailingOnly = T)

if (length(args) != 1){
    stop("Please provide the config file or --version")
}
if (args[1] == "--version"){
    message(VERSION)
    quit()
}


## -------- Valid Operations -------- #
suppressMessages(require(methods))
validops.matrix.general <- c('colSums', 'rowSums', 'apply')
validops.vector.numeric <- c('mean', 'max', 'min', 'sum', 'median','summary', 'sd')
validops.vector.general <- c('mapply', 'lapply', 'sapply')
validops.oneVal.numeric <- c('log', 'log10', 'abs', 'exp', 'floor', 'round', 'sqrt', 'ceiling')
validops.oneVal.general <- c('as.hexmode', 'as.logical', 'as.character', 'as.octmode', 'as.integer', 'as.numeric', 'toString')
validops.twoVal.numeric <- c('^', '%/%', '%%', '*', '-', '+', '/')
validops.twoVal.boolean <- c('<', '<=', '>', '>=', '==', '!=')
validops.twoVal.boolean.map <- list(lt="<", leq="<=", gt=">", geq=">=", eq="==", neq="!=")    # convenience map for Galaxy

## Random utils for sub-selecting from a data.frame
validops.dataframes <- c('[.AsIs', '[.data.frame', '[.Date', '[.difftime', '[.Dlist', '[.factor', '[.hexmode', '[.listof', '[.noquote', '[.numeric_version', '[.octmode', '[.POSIXct', '[.POSIXlt', '[.simple.list', '[.table', '[.warnings', '[<-', '[<-.data.frame', '[<-.Date', '[<-.factor', '[<-.numeric_version', '[<-.POSIXct', '[<-.POSIXlt', '[[', '[[.data.frame', '[[.Date', '[[.factor', '[[.numeric_version', '[[.POSIXct', 'Math.data.frame', 'as.matrix', c(methods(class=data.frame)))


## General operations
validops.general <- c('^', '{', '|', '|.hexmode', '|.octmode', "I", "=", "$", "<<-","<-", "{", "(", "[",":", "if", "order", "function", "print", "c", "rev", "rev.default", "print.default", "[.data.frame", "sort", "unique", "sort.default", "unique.default", "paste", "parse", "eval", "is.null", "length", "hist", "hist.default", "trimws", "return", "names", "names<-", "!", "is.na", "list", "mean", "mean.default", "sd", "assign", "data.frame", "integer", "character", "logical", "rbind", "cbind", "options")



## -------- Valid Operations: Methods -------- #
getTwoValueBooleanOperator <- function(op){
    #' Returns boolean function that compares one value with another
    #'
    #' @param op the name of the operation
    #' @return a function
    if (op %in% validops.twoVal.boolean) return(get(op))
    ## Galaxy convenience map
    trans <- names(validops.twoVal.boolean.map)
    if (op %in% trans) return(get(validops.twoVal.boolean.map[[op]]))
    stop(op, " is not a valid two value boolean operator")
}

getTwoValueNumericOperator <- function(op){
    #' Returns numeric function that operates one value upon another
    #'
    #' @param op the name of the operation
    #' @return a function
    if (op %in% validops.twoVal.numeric) return(get(op))
    stop(op, " is not a two value numeric operator")
}

getOneValueOperator <- function(op){
    #' Returns a unary operator
    #'
    #' @param op the name of the operation
    #' @return a function
    if (op %in% c(validops.oneVal.numeric, validops.oneVal.general)) return(get(op))
    stop(op, " is not a valid one value operator")
}

getVectorNumericOperator <- function(op){
    #' Returns a function that operates over vectors and returns a numeric value
    #'
    #' @param op the name of the operation
    #' @return a function
    if (op %in% validops.vector.numeric) return(get(op))
    stop(op, " is not a vector operator")
}

getMatrixOperator <- function(op){
    #' Returns a function that operates over 2D matrices and returns either a vector
    #' or a matrix.
    #'
    #' @param op the name of the operation
    #' @return a function
    if (op %in% validops.matrix.general) return(get(op))
    stop(op, " in not a matrix operator")
}

## -------- Safe Environment Class Definition -------- #
##
## Here we define our safe environment for the evaluation of general table functions
##
SafeEnv <- setClass("SafeEnv", slots= c(env = "environment", safelist = "vector" ))

## Quick primer in R S4 object 'classes':
## - initialize is the default constructor
## - setGeneric makes the function visible to R, outlines a template of args
## - setMethod binds (and redefines the fn body) to any class of the same signature.
setMethod(
    "initialize",
    signature("SafeEnv"),
    definition = function(.Object, safelist){
        ## setValidity methods are generally frowned upon since the constructor should already do this
        if (! is.vector(safelist)){
            stop("safelist is not a vector!")
        }
        .Object@env <- new.env(parent = emptyenv())
        ## Packages should only be loaded on init
        for (f in safelist){
            tryCatch({
                .Object@env[[f]] <- get(f, "package:base")
            },
            error = function(e) {
                tryCatch({
                    .Object@env[[f]] <- get(f)
                },
                error = function(e){
                    ## Warnings, not fatal
                    write(paste("Warning: Could not load", f), stdout())
                })
            })
        }
        return(.Object)
    }
)
setGeneric("addToEnv", function(.Object,name,new_obj){standardGeneric("addToEnv")})
setMethod(
    "addToEnv",
    signature("SafeEnv"),
    function(.Object, name, new_obj){
        if (exists(name, envir = .Object@env)){
            stop(name, "already exists in safe environment")
        } else {
            .Object@env[[name]] <- new_obj
        }
    }
)
setGeneric("safeget", function(object,x) standardGeneric("safeget") )
setMethod(
    "safeget",
    signature = "SafeEnv",
    definition = function(object, x){
        eval(substitute(x), env = object@env)
    }
)
##
## Initiate our safe environment
##
se <- SafeEnv(c(
    getGroupMembers("Math"),
    getGroupMembers("Math2"),
    getGroupMembers("Arith"),
    getGroupMembers("Compare"),
    getGroupMembers("Summary"),
    getGroupMembers("Ops"),
    getGroupMembers("Complex"),
    validops.dataframes,
    validops.general,
    validops.twoVal.boolean,
    validops.twoVal.numeric,
    validops.oneVal.numeric,
    validops.oneVal.general,
    validops.vector.numeric,
    validops.vector.general,
    validops.matrix.general
))




## -------- Parsing the Config File -------- #
sanityCheck <- function(file){
    #' Perform an initial sanityCheck for obvious signs of malicious intent
    #'
    #' We cannot use the 'source' function in our safe environment, since users could
    #' load libraries at will. Instead we will parse the file line by line to try to
    #' block any suspicious file or library requests.
    #'
    #' This check must be performed in the global environment before the actual
    #' evaluation of code is performed in the safe environment
    #'
    #' @param file input config file provided by Galaxy
    #' @return lines from the config file
    lines <- readLines(file, warn=F)
    for (lin in lines){
        if (
            grepl("function", lin) ||
            grepl("\\.read|read[.(]", lin) ||
            grepl("\\.write|write[.(]", lin) ||
            grepl("^file|file[.(]", lin) ||
            grepl("list.files", lin) ||
            grepl("source\\(", lin) ||
            grepl("require\\(", lin)
        ){
            stop("Illegal call: ", lin)
        }
    }
    return(lines)
}

zzz <- addToEnv(se, "string2vector", function(string){
    #' Converts a string to a valid vector
    #'
    #' e.g. "-1,2:5,2"  evaluates to c(-1,2,3,4,5,2)
    #'
    #' @param add string2vector function to safe environment
    #' @return vector of integers
    terms = unlist(strsplit(string, split=","))
    res <- lapply(terms, function(x){
        if (grepl(":",x)){
            l_r <- as.integer(c(unlist(strsplit(x,split=":"))))
            return(seq(l_r[1],l_r[2]))
        }
        return(as.integer(x))
    })
    return(c(unlist(res)))
})

##
## Add our config file to the safe environment, if it passes the initial check
##
zzz <- addToEnv(se, "conf", sanityCheck(args[1]))

##
## Evaluate the config text in the safe environment.
##
zzz <- safeget(se, eval(parse(text=conf)))


## -------- HERE BE QUIET WATERS -------- #
##
## At this point our safe environment has been deemed 'safe'
##


##
## Since file operations are limited to the global environment, we need to load the input table
## into our global environment (using the settings from the safe environment) and only then load
## it back to our safe environment.
##
reader.skip <- safeget(se, reader.skip)
customfunc.header <- NULL

logfile <- stderr()
outtable <- safeget(se, outtable)

##
## Load table into safe environment, and initialise empty output table variable
##
multiple.tables <- safeget(se, multiple.tables)

if (length(multiple.tables) > 0){
    for (i in 1:length(multiple.tables)){
        tdat <- multiple.tables[i,]
        if (!is.na(tdat$file)){
            tmp <- read.table(tdat$file, header=tdat$header, row.names=tdat$row.names, skip=tdat$skip)
            if (nrow(tmp)==1 || ncol(tmp)==1){
                tmp <- unlist(c(tmp))
            }
            else if (tdat$ismatrix) {
                tmp <- as.matrix(tmp)
            }
            zzz <- addToEnv(se, paste("table",i,sep=""), tmp)
        }
    }
} else {
    reader.file <- safeget(se, reader.file)
    reader.header <- safeget(se, reader.header)
    reader.row.col <- safeget(se, reader.row.col)
    reader.ismatrix <- safeget(se, reader.ismatrix)

    tmp <- read.table(reader.file, header = reader.header, row.names = reader.row.col, skip = reader.skip)
    if (reader.ismatrix){
        tmp <- as.matrix(tmp)
    }
    zzz <- addToEnv(se, "tab", tmp)
}
zzz <- addToEnv(se, "out", c())


## -------- Main -------- #
user.mode <- safeget(se, user.mode)

if (user.mode == "whole_table_operations"){
    zzz <- safeget(
        se,
        x = {
            eval(parse(
                text = paste(
                    "tableop <- function(){ ", fulltable.customop, " }"
                )
            ))
            out <- tableop()
        }
    )
    customfunc.header <- safeget(se, fulltable.customop)

    writeLines(c(
        "Multiple Table Mode:",
        safeget(se, fulltable.customop)
    ), logfile)

} else if (user.mode == "sort") {
    ## Sorting using column names is easy, but using column indices
    ## requires passing vectors into the order() command as a non-list or non-vector
    ##    e.g. dd[,(order(dd[1,], -dd[3,])]
    ##
    ## However, this requires knowing the number of indices beforehard, and unlist()
    ## does not appear to work here. There are examples for this scenario in the help
    ## text, but they all apply to named vectors.
    ##
    ## An easier way to is to simply pass successive single vectors in the
    ## the opposite order they are given.
    ##    e.g.  if we want column 3 in descending and 1 in ascending order:
    ##               dd <- dd[,order(dd[1,])];  dd <- dd[,order(-dd[3,])];
    zzz <- addToEnv(se, "dd", c())
    zzz <- safeget(
        se,
        x = {
            dd <- tab

            if (sort.mode == "rows"){  # sort rows by column number
                res <- sapply(rev(sort.fields), FUN=function(a){
                    if (a < 0){
                        dd <<- dd[order(-dd[,-a]),]
                    } else {
                        dd <<- dd[order(dd[,a]),]
                    }
                })
            }

            else if (sort.mode == "cols"){ # sort cols by row number
                res <- sapply(rev(sort.fields), FUN=function(a){
                    if (a < 0) {
                        ##dd <<- dd[,order(-dd[-a,])]
                        ##
                        ## This is genuinely NUTS. There is a specific function in R that
                        ## permits:
                        ##    "-dd[,-a]"     but NOT      "-dd[-a,]"
                        ##
                        ## I will waste no more time looking for this, and for now will
                        ## just make do with this horrible hack:
                        dd <<- dd[,order(apply(dd, 1:2, function(x){ x * -1 })[-a,])]
                    } else {
                        dd <<- dd[,order(dd[a,])]
                    }
                })
            }
            out <- dd
        }
    )
    writeLines(c(
        paste("Sort Mode: ",
              safeget(se, sort.mode), safeget(se, sort.fields)
              )), logfile)

} else if (user.mode == "select"){
    zzz <- safeget(
        se,
        x = {
            if (select.cols.unique){
                select.cols.wanted <- unique(sort(select.cols.wanted))
            }

            if (select.rows.unique){
                select.rows.wanted <- unique(sort(select.rows.wanted))
            }
            out <- tab[select.rows.wanted,select.cols.wanted]
        }
    )
    writeLines(c(
        paste("Select Mode:",
              safeget(se, select.cols.unique), safeget(se, select.cols.wanted),
              safeget(se, select.rows.unique), safeget(se, select.rows.wanted)
              )), logfile)



} else if (user.mode == "filtersum"){
    fop <- getTwoValueBooleanOperator(safeget(se, filtersum.op))
    zzz <- addToEnv(se, "fop", fop)
    zzz <- addToEnv(se, "vec", c())
    zzz <- safeget(
        se,
        x = {
            if (filtersum.mode == "rowSum"){
                vec <- rowSums(tab, na.rm = narm)
                out <- tab[fop(vec, filtersum.value),]
            } else if (filtersum.mode == "colSum"){
                vec <- colSums(tab, na.rm = narm)
                out <- tab[,fop(vec, filtersum.value)]
            }
        }
    )
    writeLines(c(
        paste("Filter Sum Mode: ",
              safeget(se, filtersum.mode),
              safeget(se, filtersum.op),
              safeget(se, filtersum.value),
              "- Values modified"
              )), logfile)
    writeLines(paste(safeget(se, vec)), logfile)

} else if (user.mode == "matrixapply"){
    use_custom <- safeget(se, matrixapply.custom)

    if (!use_custom){
        zzz <- addToEnv(se, "op",
                 getVectorNumericOperator(
                     safeget(se,matrixapply.op))
                 )
        zzz <- safeget(
            se,
            x = {
                out <- apply(tab, matrixapply.dimension, op)
            }
        )
    } else {
        ##
        ## Custom function
        ##
        zzz <- safeget(
            se,
            x = {
                fun_body <- paste(
                    "op <- function(vec){ ", matrixapply.custom.func, " }"
                )
                eval(parse(text = fun_body))
                out <- apply(tab, matrixapply.dimension, op)
            }
        )
        zzz <- safeget(se, matrixapply.op <- matrixapply.custom.func)
    }

    customfunc.header <- safeget(se, matrixapply.op)
    
    writeLines(c(
        safeget(se, paste(
                         "Matrix Mode: ",
                         matrixapply.op, matrixapply.dimension
                     ))
    ), logfile)
} else if (user.mode == "element"){
    zzz <- addToEnv(se, "getTwoValueBooleanOperator", getTwoValueBooleanOperator)
    zzz <- addToEnv(se, "getTwoValueNumericOperator", getTwoValueNumericOperator)
    zzz <- addToEnv(se, "getOneValueOperator", getOneValueOperator)

    zzz <- safeget(
        se,
        x = {
            ## NOTE:
            ##  Unable to find the base function that allows for simple map to happen:
            ##        tab[change_elems] <- op2(tab[change_elems], element.scale.value)
            ##
            ## For now, we will just use apply.
            ##
            ##
            ## Sub-selection of elements to modify occurs only if change_elems is not true
            ##
            op = TRUE
            change_elems = TRUE

            if (element.op != "T"){
                val <- element.value
                val <- if (is.na(as.numeric(val))) val else as.numeric(val)

                op <- getTwoValueBooleanOperator(element.op)
                change_elems <- apply(tab, 1:2, function(x){
                    op(x, val)
                })
            }

            if (element.mode == "replace") {
                val <- element.replace
                val <- if (is.na(as.numeric(val))) val else as.numeric(val)

                tab[change_elems] <- sapply(tab[change_elems], function(x){
                    val
                })
            } else if (element.mode == "modify"){
                op2 <- getOneValueOperator(element.modify.op)
                tab[change_elems] <- sapply(tab[change_elems], function(x){
                    op2(x)
                })
            }  else if (element.mode == "scale"){
                op2 <- getTwoValueNumericOperator(element.scale.op)
                tab[change_elems] <- sapply(tab[change_elems], function(x){
                    op2(x,element.scale.value)
                })
            } else if (element.mode == "custom"){
                ##
                ## Custom function
                ##
                eval(parse(
                    text = paste(
                        "elemop <- function(elem){ ", element.customop, " }"
                    )
                ))
                tab[change_elems] <- sapply(tab[change_elems], elemop)
            }
            out <- tab
        }
    )

    writeLines(c(
        safeget(se, paste("Elements Mode: ", element.mode, element.op, element.value)),
        safeget(se, paste("- Values modified:", if (is.null(element.modify.op)) "" else paste(element.modify.op, element.modify.op))),
        safeget(se, if (length(change_elems) > 1) change_elems else "all")
    ), logfile)

} else if (user.mode == "fulltable"){
    ##
    ## Entire section is a custom function
    ##
    zzz <- safeget(
        se,
        x = {
            eval(parse(
                text = paste(
                    "tableop <- function(table){ ", fulltable.customop, " }"
                )
            ))
            out <- tableop(tab)
        }
    )
    customfunc.header <- safeget(se, fulltable.customop)

    
    writeLines(c(
        "Full Table Mode:",
        safeget(se, fulltable.customop)
    ), logfile)
}


if (user.mode != "histogram"){
    ##
    ## Write out table data
    ##
    out <- safeget(se, out)

    is.vector <- !(class(out) %in% c("data.frame", "matrix"))
    is.custom <- !(is.null(customfunc.header))
    is.singtb <- (!is.vector && min(dim(out))==1)

       ## If blank func, or output table is a not vector
    if (is.custom && (is.vector || is.singtb)){
        ## First line of outtable is function name
        ## - This is a hack since list(eval(matrixapply.custom.func)=out) is not
        ##   possible. Most characters are replaced with dots due to write.table limitations
        write(paste("\t", customfunc.header, sep=""), outtable)
        write.table(out, outtable, sep='\t', quote=FALSE, col.names=FALSE, append=TRUE)
    } else {
        write.table(out, outtable, sep='\t', quote=FALSE, col.names=NA, append=FALSE)
    }
} else {
    ##
    ## Draw histograms
    ##
    zzz <- safeget(se, x = {

        colvals <- colSums(tab, na.rm = narm)
        rowvals <- rowSums(tab, na.rm = narm)
        xlabtext = "Sums"

        if (histo.mod == "ln"){
            colvals <- log(colvals)
            rowvals <- log(rowvals)
            xlabtext = "ln(Sums)"

        } else if (histo.mod == "log10"){
            colvals <- log10(colvals)
            rowvals <- log10(rowvals)
            xlabtext = "log10(Sums)"
        }

        colxlims = if (is.null(histo.custom_limits.cols.x)) c(min(colvals, na.rm=narm),max(colvals, na.rm=narm)) else histo.custom_limits.cols.x
        rowxlims = if (is.null(histo.custom_limits.rows.x)) c(min(rowvals, na.rm=narm),max(rowvals, na.rm=narm)) else histo.custom_limits.rows.x
        colylims = histo.custom_limits.cols.y
        rowylims = histo.custom_limits.rows.y
    })

    png(safeget(se,outhisto))
    par(mfrow=c(1,2))

    zzz <- safeget(se, x = {
        hist(colvals, breaks = histo.breaks.cols, xlab = xlabtext, xlim = colxlims, ylim = colylims,
             main = "Histogram of Column Sums")
        hist(rowvals, breaks = histo.breaks.rows, xlab = xlabtext, xlim = rowxlims, ylim = rowylims,
             main = "Histogram of Row Sums" )
    })
    dev.off()

    writeLines(c(
        paste("Histogram Mode: ", safeget(se,histo.mod))
    ), logfile)
}
