#!/usr/bin/env Rscript

# Author: Etienne CAMENEN
# Date: 2021
# Contact: etienne.camenen@gmail.com
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

rm(list = ls())
graphics.off()
separator <- NULL

########## Arguments ##########

# Parse the arguments from a command line launch
get_args <- function() {
    option_list <- list(
        # File parameters
        make_option(
            opt_str = c("-d", "--datasets"),
            type = "character",
            metavar = "path list",
            help = "List of comma-separated file paths corresponding to the
            blocks to be analyzed (one per block and without spaces between
            them; e.g., path/file1.txt,path/file2.txt) [required]"
        ),
        make_option(
            opt_str = c("-c", "--connection"),
            type = "character",
            metavar = "path",
            help = "Path of the file defining the connections between the blocks
            [if not used, activates the superblock mode]"
        ),
        make_option(
            opt_str = "--group",
            type = "character",
            metavar = "path",
            help = "Path of the file coloring the individuals in the ad hoc
            plot"
        ),
        make_option(
            opt_str = c("-r", "--response"),
            type = "integer",
            metavar = "integer",
            help = "Position of the response file for the supervised mode within
            the block path list [actives the supervised mode]"
        ),
        make_option(
            opt_str = "--names",
            type = "character",
            metavar = "character list",
            help = "List of comma-separated block names to rename them (one per
            block; without spaces between them) [default: the block file names]"
        ),
        make_option(
            opt_str = c("-H", "--header"),
            type = "logical",
            action = "store_false",
            help = "DO NOT consider the first row as the column header"
        ),
        make_option(
            opt_str = "--separator",
            type = "integer",
            metavar = "integer",
            default = opt[1],
            help = "Character used to separate columns (1: tabulation,
            2: semicolon, 3: comma) [default: %default]"
        ),
        # Analysis parameter
        make_option(
            opt_str = "--type",
            type = "character",
            metavar = "character",
            default = opt[2],
            help = "Type of analysis [default: %default] (among: rgcca, sgcca,
            pca, spca, pls, spls, cca, ifa, ra, gcca, maxvar, maxvar-b,
            maxvar-a, mcoa,cpca-1, cpca-2, cpca-4, hpca, maxbet-b, maxbet,
            maxdiff-b, maxdiff, maxvar-a, sabscor, ssqcor, ssqcov-1, ssqcov-2,
            ssqcov, sumcor, sumcov-1, sumcov-2, sumcov, sabscov, sabscov-1,
            sabscov-2)"
        ),
        make_option(
            opt_str = "--ncomp",
            type = "character",
            metavar = "integer list",
            default = opt[3],
            help = "Number of components in the analysis for each block
            [default: %default]. The number should be higher than 1 and lower
            than the minimum number of variables among the blocks. It can be a
            single values or a comma-separated list (e.g 2,2,3,2)."
        ),
        make_option(
            opt_str = "--penalty",
            type = "character",
            metavar = "float list",
            default = opt[4],
            help = "For RGCCA, a regularization parameter for each block (i.e., tau)
            [default: %default]. Tau varies from 0 (maximizing the correlation)
            to 1 (maximizing the covariance). For SGCCA, tau is automatically
            set to 1 and shrinkage parameter can be defined instead for
            automatic variable selection, varying from the square root of the
            variable number (the fewest selected variables) to 1 (all the
            variables are included). It can be a single value or a
            comma-separated list (e.g. 0,1,0.75,1)."
        ),
        make_option(
            opt_str = "--scheme",
            type = "integer",
            metavar = "integer",
            default = opt[5],
            help = "Link (i.e. scheme) function for covariance maximization
            (1: x, 2: x^2, 3: |x|, 4: x^4) [default: %default]. Onnly, the x
            function ('horst scheme') penalizes structural negative correlation.
            The x^2 function ('factorial scheme') discriminates more strongly
            the blocks than the |x| ('centroid scheme') one."
        ),
        make_option(
            opt_str = "--scale",
            type = "logical",
            action = "store_false",
            help = "DO NOT scale the blocks (i.e., a data centering step is
            always performed). Otherwise, each block is normalised and divided
            by the squareroot of its number of variables."
        ),
        make_option(
            opt_str = "--superblock",
            type = "logical",
            action = "store_false",
            help = "DO NOT use a superblock (i.e. a concatenation of all the
            blocks to visualize them all together in a consensus space). In
            this case, all blocks are assumed to be connected or a connection
            file could be used."
        ),
        # Graphical parameters
        make_option(
            opt_str = "--text",
            type = "logical",
            action = "store_false",
            help = "DO NOT display the name of the points instead of shapes when
            plotting"
        ),
        make_option(
            opt_str = "--block",
            type = "integer",
            metavar = "integer",
            default = opt[6],
            help = "Position in the path list of the plotted block (0: the
            superblock or, if not activated, the last one, 1: the fist one,
            2: the 2nd, etc.)[default: the last one]"
        ),
        make_option(
            opt_str = "--block_y",
            type = "integer",
            metavar = "integer",
            help = "Position in the path list of the plotted block for the
            Y-axis in the individual plot (0: the superblock or, if not
            activated, the last one, 1: the fist one, 2: the 2nd, etc.)
            [default: the last one]"
        ),
        make_option(
            opt_str = "--compx",
            type = "integer",
            metavar = "integer",
            default = opt[7],
            help = "Component used in the X-axis for biplots and the only
            component used for histograms [default: %default] (should not be
            higher than the number of components of the analysis)"
        ),
        make_option(
            opt_str = "--compy",
            type = "integer",
            metavar = "integer",
            default = opt[8],
            help = "Component used in the Y-axis for biplots
            [default: %default] (should not be higher than the number of
            components of the analysis)"
        ),
        make_option(
            opt_str = "--nmark",
            type = "integer",
            metavar = "integer",
            default = opt[9],
            help = "Number maximum of top variables in ad hoc plot
            [default: %default]"
        ),
        # output parameters
        make_option(
            opt_str = "--o1",
            type = "character",
            metavar = "path",
            default = opt[10],
            help = "Path for the individual plot [default: %default]"
        ),
        make_option(
            opt_str = "--o2",
            type = "character",
            metavar = "path",
            default = opt[11],
            help = "Path for the variable plot [default: %default]"
        ),
        make_option(
            opt_str = "--o3",
            type = "character",
            metavar = "path",
            default = opt[12],
            help = "Path for the top variables plot [default: %default]"
        ),
        make_option(
            opt_str = "--o4",
            type = "character",
            metavar = "path",
            default = opt[13],
            help = "Path for the explained variance plot [default: %default]"
        ),
        make_option(
            opt_str = "--o5",
            type = "character",
            metavar = "path",
            default = opt[14],
            help = "Path for the design plot [default: %default]"
        ),
        make_option(
            opt_str = "--o6",
            type = "character",
            metavar = "path",
            default = opt[15],
            help = "Path for the individual table [default: %default]"
        ),
        make_option(
            opt_str = "--o7",
            type = "character",
            metavar = "path",
            default = opt[16],
            help = "Path for the variable table [default: %default]"
        ),
        make_option(
            opt_str = "--o8",
            type = "character",
            metavar = "path",
            default = opt[17],
            help = "Path for the analysis results in RData [default: %default]"
        )
    )
    return(optparse::OptionParser(option_list = option_list))
}

check_arg <- function(opt) {
    # Check the validity of the arguments opt : an optionParser object

    if (is.null(opt$datasets)) {
        stop_rgcca(paste0("datasets is required."), exit_code = 121)
    }

    if (is.null(opt$scheme)) {
        opt$scheme <- "factorial"
    } else if (!opt$scheme %in% seq(4)) {
        stop_rgcca(
            paste0(
                "scheme should be comprise between 1 and 4 [by default: 2], not be equal to ",
                opt$scheme,
                "."
            ),
            exit_code = 122
        )
    } else {
        schemes <- c("horst", "factorial", "centroid")
        if (opt$scheme == 4) {
            opt$scheme <- function(x) x^4
        } else {
            opt$scheme <- schemes[opt$scheme]
        }
    }

    if (!opt$separator %in% seq(3)) {
        stop_rgcca(
            paste0(
                "separator should be comprise between 1 and 3 (1: Tabulation, 2: Semicolon, 3: Comma) [by default: 2], not be equal to ",
                opt$separator,
                "."
            ),
            exit_code = 123
        )
    } else {
        separators <- c("\t", ";", ",")
        opt$separator <- separators[opt$separator]
    }

    nmark <- NULL
    RGCCA:::check_integer("nmark", opt$nmark, min = 2)

    for (x in c("ncomp", "penalty")) {
        opt[[x]] <- char_to_list(opt[[x]])
    }

    return(opt)
}

post_check_arg <- function(opt, rgcca) {
    # Check the validity of the arguments after loading the blocks opt : an
    # optionParser object blocks : a list of matrix
    blocks <- NULL
    for (x in c("block", "block_y")) {
        if (!is.null(opt[[x]])) {
            if (opt[[x]] == 0) {
                opt[[x]] <- length(rgcca$call$blocks)
            }
            opt[[x]] <- RGCCA:::check_blockx(x, opt[[x]], rgcca$call$blocks)
        }
    }

    if (any(opt$ncomp == 1)) {
        opt$compy <- 1
    }

    for (x in c("compx", "compy")) {
        opt[[x]] <- check_compx(x, opt[[x]], rgcca$call$ncomp, opt$block)
    }

    return(opt)
}

########## Main ##########

# Get arguments : R packaging install, need an opt variable with associated
# arguments
opt <- list(
    separator = 1,
    type = "rgcca",
    ncomp = 2,
    penalty = 1,
    scheme = 2,
    block = 0,
    compx = 1,
    compy = 2,
    nmark = 100,
    o1 = "individuals.pdf",
    o2 = "corcircle.pdf",
    o3 = "top_variables.pdf",
    o4 = "ave.pdf",
    o5 = "design.pdf",
    o6 = "individuals.tsv",
    o7 = "variables.tsv",
    o8 = "rgcca_result.RData",
    datasets = paste0("inst/extdata/",
        c("agriculture", "industry", "politic"),
        ".tsv",
        collapse = ","
    )
)

# Load functions
all_funcs <- unclass(lsf.str(envir = asNamespace("RGCCA"), all = TRUE))
for (i in all_funcs) {
    eval(parse(text = paste0(i, "<-RGCCA:::", i)))
}

load_libraries(c("ggplot2", "optparse", "scales", "igraph", "MASS", "Deriv"))
try(load_libraries("ggrepel"), silent = TRUE)

tryCatch(
    opt <- check_arg(optparse::parse_args(get_args())),
    error = function(e) {
        if (length(grep("nextArg", e[[1]])) != 1) {
            stop_rgcca(e[[1]], exit_code = 140)
        }
    }, warning = function(w) {
        stop_rgcca(w[[1]], exit_code = 141)
    }
)

# Set missing parameters by default
opt$header <- !("header" %in% names(opt))
opt$superblock <- !("superblock" %in% names(opt))
opt$scale <- !("scale" %in% names(opt))
opt$text <- !("text" %in% names(opt))
cex_lab <- 20
cex_main <- 25
cex_point <- 3
cex_sub <- 20
cex_axis <- 10
cex <- 1.25

status <- 0
tryCatch(
    {
        blocks <- load_blocks(opt$datasets, opt$names, opt$separator)
        group <- load_response(blocks, opt$group, opt$separator, opt$header)
        connection <- load_connection(file = opt$connection, separator = opt$separator)

        func <- quote(
            rgcca(
                blocks = blocks,
                connection = connection,
                response = opt$response,
                superblock = opt$superblock,
                ncomp = opt$ncomp,
                scheme = opt$scheme,
                scale = opt$scale,
                method = opt$type
            )
        )
        if (tolower(opt$type) %in% c("sgcca", "spca", "spls")) {
            func[["sparsity"]] <- opt$penalty
        } else {
            func[["tau"]] <- opt$penalty
        }

        rgcca_out <- eval(as.call(func))

        opt <- post_check_arg(opt, rgcca_out)

        ########## Plot ##########

        if (rgcca_out$call$ncomp[opt$block] == 1 && is.null(opt$block_y)) {
            warning("With a number of component of 1, a second block should be chosen to perform an individual plot")
        } else {
            (
                individual_plot <- plot_ind(
                    rgcca_out,
                    group,
                    opt$compx,
                    opt$compy,
                    opt$block,
                    opt$text,
                    opt$block_y,
                    "Response",
                    cex_lab = cex_lab,
                    cex_point = cex_point,
                    cex_main = cex_main,
                    cex = cex
                )
            )
            save_plot(opt$o1, individual_plot)
        }

        if (rgcca_out$call$ncomp[opt$block] > 1) {
            (
                corcircle <- plot_var_2D(
                    rgcca_out,
                    opt$compx,
                    opt$compy,
                    opt$block,
                    opt$text,
                    n_mark = opt$nmark,
                    cex_lab = cex_lab,
                    cex_point = cex_point,
                    cex_main = cex_main,
                    cex = cex
                )
            )
            save_plot(opt$o2, corcircle)
        }

        top_variables <- plot_var_1D(
            rgcca_out,
            opt$compx,
            opt$nmark,
            opt$block,
            type = "loadings",
            title = paste0("Variable correlations", ": ", names(rgcca_out$call$blocks)[opt$block], " with "),
            cex_sub = cex_sub,
            cex_main = cex_main,
            cex_axis = cex_axis,
            cex = cex
        )
        save_plot(opt$o3, top_variables)

        # Average Variance Explained
        (ave <- plot_ave(
            rgcca_out,
            cex_main = cex_main,
            cex_sub = cex_sub,
            cex_axis = cex_axis,
            cex = cex
        ))
        save_plot(opt$o4, ave)

        # Creates design scheme
        design <- function() {
            plot_network(
                rgcca_out,
                cex_main = cex_main,
                cex_point = cex_point,
                cex = cex
            )
        }
        save_plot(opt$o5, design)

        save_ind(rgcca_out, opt$o6)
        save_var(rgcca_out, opt$o7)
        save(rgcca_out, file = opt$o8)
    },
    error = function(e) {
        if (class(e)[1] %in% c("simpleError", "error", "condition")) {
            status <<- 1
        } else {
            status <<- class(e)[1]
        }
        msg <- "The design matrix C"
        if (grepl(msg, e$message)) {
            e$message <- gsub(msg, "The connection file", e$message)
        }
        message(e$message)
    }
)
quit(status = status)
