# libraries
suppressMessages(library(ggfortify))
suppressMessages(library(GetoptLong))
suppressMessages(library(scales))

# defaults VARs
png_res_width <- 1600
png_res_height <- 1400
svg_res_width <- 16
svg_res_height <- 14
output2 <- character(0)

# static VARs
groups <- TRUE
circles_boolean <- TRUE
group_name <- ""
group_cols <- ""
group_feature <- ""
group_colors <- ""
group_colors_vector <- ""

# args <- commandArgs(trailingOnly = TRUE)

# test VARS
# 
# input <- "~/projects/qPCR_plots/test-data/example1.txt"
# header <- TRUE
# rowname_index <- 1
# horizontal <- TRUE
# title <- "Test"
# transform <- "none"
# background <- "Default"
# scaling <- "Automatic"
# legend <- TRUE
# group_type <- "define_groups"
# group_name <- c("ff", "ggg", "zzz", "jj")
# group_cols <- c("1,2,3,4,5", "6,7,8,9,10", "11,12,13,14,15", "16,17,18,19,20")
# group_colors <- c("none", "#ffc000", "none", "#76923c")
# plot_param <- "label"
# circles <- "convex"
# xaxismin <- 0
# xaxismax <- 1
# yaxismin <- 2 
# yaxismax <- 3
# log_modi <- 0
# png_res_height <- 1400
# png_res_width <- 1600
# png_res <- 200
# output_type <- FALSE
# svg_res_height <- 14
# svg_res_width <- 16
# output1 <- "out.png"
# output2 <- "None"

# inputs
GetoptLong(
  "input=s", "Input file, string (tabular), mandatory option",
  "header!", "header, boolean, mandatory option",
  "rowname_index=i", "ggplots title, integer, mandatory option",
  "horizontal!", "horizonal sample names, boolean, mandatory option",
  "title=s", "ggplots title, string, mandatory option",
  "transform=s", "transform data to log2 log2+1 log10 log10+1, string, mandatory option",
  "scaling=s", "scaling of axis, string, mandatory option",
  "background=s", "background of plot, string, mandatory option",
  "legend!", "legend of dataset, boolean, mandatory option",
  "xaxismin=f", "xaxismin, float, mandatory option",
  "xaxismax=f", "xaxismax, float, mandatory option",
  "yaxismin=f", "yaxismin, float, mandatory option",
  "yaxismax=f", "yaxismax, float, mandatory option",
  "group_type=s", "Group type, string, mandatory option",
  "group_name=s@", "Group name, string, optional",
  "group_cols=s@", "Group colums, string, optional",
  "group_colors=s@", "Group colors, string, optional",
  "circles=s", "Elyptic parameter - (no, convex, t, norm, euclid), string, mandatory option",
  "plot_param=s", "Plot parameter - comma separated string, string, mandatory option",
  "png_res_width=i", "PNG width , integer, optional (default 1600px)",
  "png_res_height=i", "PNG height , integer, optional (default 1400px)",
  "png_res=i", "PNG resolution , integer, optional (default 200)",
  "output_type!", "output type (SVG yes, no?), boolean, mandatory option",
  "svg_res_width=i", "SVG width , integer, optional (default 16in)",
  "svg_res_height=i", "SVG height , integer, optional (default 14in)",
  "output1=s", "Output1 PNG file path, string, mandatory option",
  "output2=s", "Output2 SVG file path, string, optional",
  "verbose",  "print messages",
  head = 'Rscript ggplot_pca.R ',
  foot = 'Please contact jochen.bick@usys.ethz.ch for comments'
)

# split plotoptions
plot_param_list <- strsplit(x = plot_param, split = ",")
plot_options <- c("shape", "label") %in% plot_param_list[[1]]

#cat("\ngroup_colors: ", group_colors)

# read table with or with out header or row_names
if(rowname_index > 0){
  df <- read.table(input, header = header, row.names = rowname_index, sep = "\t")
}else{
  df <- read.table(input, header = header, sep = "\t")
}

# check if indices are out of range
num_cols <- length(names(df))
check_col_indices <- as.integer(unlist(strsplit(group_cols, split = ",")))
if(any(check_col_indices>num_cols)){ 
  stop("Error: column indices for grouping are out of range! Check help!")
}

# check if table has only numbers
if(any(!sapply(df, is.numeric))){ 
  stop("Error: table contains not only numbers!")
}

# check if group_names are unique
if(length(unique(c(group_name, "no_group"))) != length(c(group_name, "no_group"))){ 
  stop("Error: group_names must be unique: ", paste(group_name, "no_group", collapse = ","), " is not unique!")
}
# stopifnot(any(!sapply(d1, is.numeric)))

# prepare group_features for grouping of samples accouring to orientation
if(horizontal){
  num_cols <- length(names(df))
  group_feature <- rep("no_group", num_cols)
}else{
  num_rows <- nrow(df)
  group_feature <- rep("no_group", num_rows)
}

# # assign colors
# num_cols <- 10
# group_name <- c("tt", "cc")
# group_cols <- c("2,3,4,5", "1")
# group_colors <- c("#548dd4", "none")

default_ggplot_colors <- ""
default_ggplot_colors_autoplot <- ""
# split group elements and assign indexes
if(group_type == "define_groups"){
  # set up colors
  color_names <- c(group_name, "no_group")
  cat("\ncolor_names: ", color_names)
  default_ggplot_colors <- hue_pal()(length(color_names))
  names(default_ggplot_colors) <- color_names
  cat("\ndefault_ggplot_colors: ", default_ggplot_colors)
  names(group_colors) <- group_name
  group_colors <- group_colors[group_colors != "none"]
  default_ggplot_colors[names(group_colors)] <- group_colors
  cat("\ndefault_ggplot_colors: ", default_ggplot_colors)
  group_string <- lapply(seq_along(group_name), function(k){
    gname <- group_name[k]
    gindex <- as.integer(strsplit(group_cols[k], split = ",")[[1]])
    gnames <- rep(gname, length(gindex)) 
    names(gindex) <- gnames
    gindex
  })
  group_string <- do.call(c, group_string)
  # if(row_names > 0){
  #   group_string <- group_string - 1
  # }
  group_feature[group_string] <- names(group_string)
  # subset colors on groups if and check if there is "no_group" present
  default_ggplot_colors <- default_ggplot_colors[unique(group_feature)]
  default_ggplot_colors_autoplot <- default_ggplot_colors[group_feature]
  cat("\ndefault_ggplot_colors: ", default_ggplot_colors)
}




# legend <- "yes"

# Print options to see what is going on
# cat("\n input: ", input)
# cat("\n title: ", title)
# cat("\n group_name: ", group_name)
# cat("\n group_cols: ", group_cols)
# cat("\n group_feature: ", group_feature)
# cat("\n transform: ", transform)

# cat("\n legend: ", legend)

# Show/hide legend
if(legend){
  gg_legend = theme(legend.position="right")
} else {
  gg_legend = theme(legend.position="none")
  cat("\n no legend")
}

# Choose between automatically scaled x and y axis or user defined
if(scaling == "Automatic"){
  gg_scalex = NULL
  gg_scaley = NULL
} else {
  gg_scalex = xlim(xaxismin, xaxismax)
  gg_scaley = ylim(yaxismin, yaxismax)
  cat("\n xaxismin: ", xaxismin)
  cat("\n xaxismax: ", xaxismax)
  cat("\n yaxismin: ", yaxismin)
  cat("\n yaxismax: ", yaxismax)
}

# Choose theme for plot
if(background == "bw"){
  gg_theme = theme_bw()
} else {
  gg_theme = NULL
}

#Choose dimensions of output pdf
# if(dimentions == "Default"){
#   gg_width = 7
#   gg_height = 7
# } else {
#   gg_width =  woutputdim
#   gg_height =  houtputdim 
# }

# transpose dataset for correct PCA alignment

# transpose data.frame for plotting if sample names are horizontal
if(horizontal){
  df <- as.data.frame(t(df))
}else{
  # nothing
}

plot_df <- df

# set group column if wanted
if(group_type %in% "no_groups"){
  groups <- FALSE
}else{  
  plot_df$group <- group_feature
  group_name <- "group"
}

# else{
#   group_row <- as.integer(group_feature)
#   if(row_names > 0){
#     group_name <- names(plot_df)[group_row-1]
#     df <- df[-(group_row-1)]
#   }else{
#     group_name <- names(plot_df)[group_row]
#     df <- df[-group_row]
#   }
# }

# set boolean elipes value to plot circle options
if(circles %in% "no"){
  circles_boolean <- FALSE
}

plot_mat <- df

# tranform dataset
if(transform == "log2"){
  plot_mat <- log2(plot_mat)
  cat("\n ", transform, " transformed")
}else if(transform == "log2plus1"){
  plot_mat <- log2(plot_mat+1)
  cat("\n ", transform, " transformed")
}else if(transform == "log10"){
  plot_mat <- log10(plot_mat)
  cat("\n ", transform, " transformed")
}else if(transform == "log10plus1"){
  plot_mat <- log10(plot_mat+1)
  cat("\n ", transform, " transformed")
}else{
  plot_mat <- plot_mat
}


# #axis text(tick) custization
# if(axistextcust == "Default"){
#   gg_axistext = theme(axis.text = element_text(color = NULL, size = NULL, face = NULL))
# } else {
#   gg_axistext = theme(axis.text = element_text(color = axistextcolor, size = axistextsize,
#                                                face = axistextface))
# }


# #plot title custimization
# if(plottitlecust == "Default"){
#   gg_plottitle = theme(plot.title = element_text(color = NULL, size = NULL, face = NULL))
# } else {
#   gg_plottitle = theme(plot.title = element_text(color = plottitlecolor, size = plottitlesize,
#                                                  face = plottitleface))
# }

# #grid line customization
# if(gridlinecust == "Default"){
#   gg_gridline = NULL
# } else if(gridlinecust == "hidemajor"){
#   gg_gridline = theme(panel.grid.major = element_blank())
# } else if(gridlinecust == "hideminor"){
#   gg_gridline = theme(panel.grid.minor = element_blank())
# } else if(gridlinecust == "hideboth"){
#   gg_gridline = theme(panel.grid.minor = element_blank(),
#   panel.grid.major = element_blank())
# } else {
# }

# plot with or without groups and set more plotting options using autoplot
if(groups){
  if(circles_boolean){
    plot_out <- autoplot(prcomp(plot_mat), data = plot_df, colour = group_name,
                         frame = T, frame.type=circles, shape = plot_options[1],
                         label = plot_options[2])
  }else{
    plot_out <- autoplot(prcomp(plot_mat), data = plot_df, colour = group_name,
                         frame = F, shape = plot_options[1], label = plot_options[2])
  }
}else{
  if(!circles_boolean){
    plot_out <- autoplot(prcomp(plot_mat), data = plot_df, frame = F,
                         shape = plot_options[1], label = plot_options[2])
  }else{
    plot_out <- autoplot(prcomp(plot_mat), data = plot_df, frame = F,
                         shape = plot_options[1], label = plot_options[2])
  }
}


# add advanced plotting options for final plot
plot_out <- plot_out +
  scale_color_manual(values=default_ggplot_colors) +
  scale_fill_manual(values=default_ggplot_colors) +
  gg_theme +
  gg_legend +
  #gg_scalex + gg_scaley +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5)) 

  # theme(plot.title=element_text(family="xkcd-Regular"), text=element_text(family="xkcd-Regular"),
  #                                               axis.text.x=element_text(family="xkcd-Regular", colour="black", size = 20),
  #                                               axis.text.y=element_text(family="xkcd-Regular", colour="black", size = 10))


# catch different output formats
png(filename = output1, width = png_res_width, height = png_res_height, res = png_res)
  plot_out
q <- dev.off()

svg(filename = output2, width = svg_res_width, height = svg_res_height)
 plot_out
q <- dev.off()

