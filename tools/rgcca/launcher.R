# Author: Etienne CAMENEN
# Date: 2020
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

rm(list = ls())
graphics.off()
separator <- nmark <- blocks <- NULL

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
            opt_str = c("-w", "--directory"),
            type = "character",
            metavar = "path",
            default = opt[1],
            help = "Path of the root folder containing the R/ (e.g. for Galaxy)
            [default: the current one]"
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
            default = 1,
            help = "Character used to separate columns (1: tabulation,
            2: semicolon, 3: comma) [default: %default]"
        ),
        # Analysis parameter
        make_option(
            opt_str = "--type",
            type = "character",
            metavar = "character",
            default = opt[3],
            help = "Type of analysis [default: %default] (among: rgcca, pca,
            cca, gcca, cpca-w, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff,
            maxvar-a, maxvar-b, maxvar, niles, r-maxvar, rcon-pca, ridge-gca,
            sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sum-pca,
            sumcor, sumcov-1, sumcov-2, sumcov)"
        ),
        make_option(
            opt_str = "--ncomp",
            type = "character",
            metavar = "integer list",
            default = opt[4],
            help = "Number of components in the analysis for each block
            [default: %default]. The number should be higher than 1 and lower
            than the minimum number of variables among the blocks. It can be a
            single values or a comma-separated list (e.g 2,2,3,2)."
        ),
        make_option(
            opt_str = "--penalty",
            type = "character",
            metavar = "float list",
            default = opt[5],
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
            default = 2,
            help = "Link (i.e. scheme) function for covariance maximization
            (1: x, 2: x^2, 3: |x|, 4: x^4) [default: %default]. Only, the x
            function penalizes structural negative correlation. The x^4
            function discriminates more strongly the blocks than the x^2 one."
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
            default = opt[8],
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
            default = opt[9],
            help = "Component used in the X-axis for biplots and the only
            component used for histograms [default: %default] (should not be
            higher than the number of components of the analysis)"
        ),
        make_option(
            opt_str = "--compy",
            type = "integer",
            metavar = "integer",
            default = opt[10],
            help = "Component used in the Y-axis for biplots
            [default: %default] (should not be higher than the number of
            components of the analysis)"
        ),
        make_option(
            opt_str = "--nmark",
            type = "integer",
            metavar = "integer",
            default = opt[11],
            help = "Number maximum of top variables in ad hoc plot
            [default: %default]"
        ),
        # output parameters
        make_option(
            opt_str = "--o1",
            type = "character",
            metavar = "path",
            default = opt[12],
            help = "Path for the variable plot [default: %default]"
        ),
        make_option(
            opt_str = "--o2",
            type = "character",
            metavar = "path",
            default = opt[13],
            help = "Path for the individual plot [default: %default]"
        ),
        make_option(
            opt_str = "--o3",
            type = "character",
            metavar = "path",
            default = opt[14],
            help = "Path for the top variables plot [default: %default]"
        ),
        make_option(
            opt_str = "--o4",
            type = "character",
            metavar = "path",
            default = opt[15],
            help = "Path for the explained variance plot [default: %default]"
        ),
        make_option(
            opt_str = "--o5",
            type = "character",
            metavar = "path",
            default = opt[16],
            help = "Path for the design plot [default: %default]"
        ),
        make_option(
            opt_str = "--o6",
            type = "character",
            metavar = "path",
            default = opt[17],
            help = "Path for the individual table [default: %default]"
        ),
        make_option(
            opt_str = "--o7",
            type = "character",
            metavar = "path",
            default = opt[18],
            help = "Path for the variable table [default: %default]"
        ),
        make_option(
            opt_str = "--o8",
            type = "character",
            metavar = "path",
            default = opt[19],
            help = "Path for the analysis results in RData [default: %default]"
        )
    )
    return(optparse::OptionParser(option_list = option_list))
}

char_to_list <- function(x) {
    strsplit(gsub(" ", "", as.character(x)), ",")[[1]]
}

check_arg <- function(opt) {
    # Check the validity of the arguments opt : an optionParser object

    if (is.null(opt$datasets))
        stop(paste0("datasets is required."), exit_code = 121)

    if (is.null(opt$scheme))
        opt$scheme <- "factorial"
    else if (!opt$scheme %in% seq(4)) {
        stop(
            paste0(
                "scheme should be comprise between 1 and 4 [by default: 2], not be equal to ",
                opt$scheme,
                "."
            ),
            exit_code = 122
        )
    } else {
        schemes <- c("horst", "factorial", "centroid")
        if (opt$scheme == 4)
            opt$scheme <- function(x) x ^ 4
        else
            opt$scheme <- schemes[opt$scheme]
    }

    if (!opt$separator %in% seq(3)) {
        stop(
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

    check_integer("nmark", opt$nmark, min = 2)

    for (x in c("ncomp", "penalty"))
        opt[[x]] <- char_to_list(opt[[x]])

    return(opt)
}

post_check_arg <- function(opt, rgcca) {
# Check the validity of the arguments after loading the blocks opt : an
# optionParser object blocks : a list of matrix

    for (x in c("block", "block_y")) {
        if (!is.null(opt[[x]])) {
            if (opt[[x]] == 0)
                opt[[x]] <- length(rgcca$call$blocks)
            opt[[x]] <- RGCCA:::check_blockx(x, opt[[x]], rgcca$call$blocks)
        }
    }

    if (opt$ncomp == 1)
        opt$compy <- 1

    for (x in c("compx", "compy"))
        opt[[x]] <- check_compx(x, opt[[x]], rgcca$call$ncomp, opt$block)

    return(opt)
}

check_integer <- function(x, y = x, type = "scalar", float = FALSE, min = 1) {

    if (is.null(y))
        y <- x

    if (type %in% c("matrix", "data.frame"))
        y_temp <- y

    y <- suppressWarnings(as.double(as.matrix(y)))

    if (any(is.na(y)))
        stop(paste(x, "should not be NA."))

    if (!is(y, "numeric"))
        stop(paste(x, "should be numeric."))

    if (type == "scalar" && length(y) != 1)
        stop(paste(x, "should be of length 1."))

    if (!float)
        y <- as.integer(y)

    if (all(y < min))
        stop(paste0(x, " should be higher than or equal to ", min, "."))

    if (type %in% c("matrix", "data.frame"))
        y <- matrix(
            y,
            dim(y_temp)[1],
            dim(y_temp)[2],
            dimnames = dimnames(y_temp)
        )

    if (type == "data.frame")
        as.data.frame(y)

    return(y)
}

load_libraries <- function(librairies) {
    for (l in librairies) {
        if (!(l %in% installed.packages()[, "Package"]))
            utils::install.packages(l, repos = "cran.us.r-project.org")
        suppressPackageStartupMessages(
            library(
                l,
                character.only = TRUE,
                warn.conflicts = FALSE,
                quietly = TRUE
        ))
    }
}

########## Main ##########

# Get arguments : R packaging install, need an opt variable with associated
# arguments
opt <- list(
    directory = ".",
    separator = "\t",
    type = "rgcca",
    ncomp = 2,
    penalty = 1,
    scheme = "factorial",
    init = 1,
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
        collapse = ",")
)

load_libraries(c("ggplot2", "optparse", "scales", "igraph", "MASS", "rlang", "Deriv"))
try(load_libraries("ggrepel"), silent = TRUE)

tryCatch(
    opt <- check_arg(parse_args(get_args())),
    error = function(e) {
        if (length(grep("nextArg", e[[1]])) != 1)
            stop(e[[1]], exit_code = 140)
    }, warning = function(w)
        stop(w[[1]], exit_code = 141)
)

# Load functions
setwd(opt$directory)

all_funcs <- unclass(lsf.str(envir = asNamespace("RGCCA"), all = T))
for (i in all_funcs)
    eval(parse(text = paste0(i, "<-RGCCA:::", i)))

# Set missing parameters by default
opt$header <- !("header" %in% names(opt))
opt$superblock <- !("superblock" %in% names(opt))
opt$scale <- !("scale" %in% names(opt))
opt$text <- !("text" %in% names(opt))

status <- 0
tryCatch({

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
            type = opt$type
        )
    )
    if (tolower(opt$type) %in% c("sgcca", "spca", "spls")) {
        func[["sparsity"]] <- opt$penalty
    }else {
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
                get_filename(opt$group)
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
                n_mark = opt$nmark
            )
        )
        save_plot(opt$o2, corcircle)
    }

    top_variables <- plot_var_1D(
            rgcca_out,
            opt$compx,
            opt$nmark,
            opt$block,
            type = "cor"
        )
    save_plot(opt$o3, top_variables)

    # Average Variance Explained
    (ave <- plot_ave(rgcca_out))
    save_plot(opt$o4, ave)

    # Creates design scheme
    design <- function() plot_network(rgcca_out)
    save_plot(opt$o5, design)

    save_ind(rgcca_out, opt$compx, opt$compy, opt$o6)
    save_var(rgcca_out, opt$compx, opt$compy, opt$o7)
    save(rgcca_out, file = opt$o8)

    }, error = function(e){
        if (class(e)[1] %in% c("simpleError", "error", "condition" ))
            status <<- 1
        else
            status <<- class(e)[1]
        message(e$message)
})
quit(status = status)
