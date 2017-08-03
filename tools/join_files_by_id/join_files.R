# libraries
library(data.table)
library(parallel)

# inputs
args <- commandArgs()

input <- gsub("--in=", "", args[grepl("--in=", args)])
header <- as.integer(gsub("--he=", "", args[grepl("--he=", args)]))
join_col <- gsub("--jc=", "", args[grepl("--jc=", args)])
separator <- gsub("--sep=", "", args[grepl("--sep=", args)])
null_char <- gsub("--nc=", "", args[grepl("--nc=", args)])
output <- gsub("--out=", "", args[grepl("--out=", args)])

# test VARS
# input <- list("test-data/df1.txt", "test-data/df2.txt", "test-data/df3.txt")
# 
# input <- list("test-data/df_big_1.txt", "test-data/df_big_2.txt", "test-data/df_big_3.txt",
#               "test-data/df_big_4.txt", "test-data/df_big_5.txt", "test-data/df_big_6.txt",
#               "test-data/df_big_7.txt", "test-data/df_big_8.txt", "test-data/df_big_9.txt",
#               "test-data/df_big_10.txt")
# header <- 1
# join_col <- 1
# separator <- "ta"
# null_char <- 0
# output <- "out"

if(header > 0){
  header <- TRUE
}else{
  header <- FALSE
}
join_col <- as.integer(join_col)

# read files into list
df_list <- lapply(input, function(x){as.data.frame(fread(x))})


#### fix the ids name for all read in tables
df_list <- lapply(df_list, function(x){
  names_x <- names(x)
  names_x[names_x == "ids"] <- "id" # to join correctly
  names_x[join_col] <- "ids"
  names(x) <- names_x
  x
})

# generate unique ids string
df0 <- lapply(df_list, function(x){
  x[join_col]
})

df0 <- data.frame(ids=unique(do.call(rbind, df0)))
df_list <- append(df0, df_list)
df_list[[1]] <- data.frame(ids=df_list[[1]]) 



ids <- df_list[[1]]
ids <- data.frame(ids = ids[order(ids$ids), "ids" ])
merged_df <- mclapply(2:length(df_list), function(x){
  merged_sub <- merge(x = ids, y = df_list[[x]], by = "ids", all.x = T, sort = F)
  merged_sub <- merged_sub[order(merged_sub$ids), ]
  merged_sub[-1]
}, mc.cores=4)

df <- cbind(ids, do.call(cbind, merged_df))

# benchmarking
# library(microbenchmark)
# microbenchmark(
#   df1 <- lapply(2:length(df_list), function(x){
#     merge(df_list[[1]], df_list[[x]], by = "ids", all.x = T)[-1]
#   }),
#   merged_df <- mclapply(2:length(df_list), function(x){
#     merge(x = ids, y = df_list[[x]], by = "ids", all.x = T, sort = F)[-1]}, mc.cores=7)
#   ,times = 4
# )

# change null_char
df[is.na(df)] <- null_char
# separator <- "ta"
# print(separator)
delim <- list(ta = "\t", do = ".", co = ",", un = "_", da = "-", sp = " ")
# print(separator)
separator <- delim[[separator]]
# write data.frame to file
write.table(df, file = output, sep = separator, row.names = F, quote = F, col.names = header)
