# Function for the data input
load_data <- function(file_name, file_extension) {
    if (file_extension == "csv") {
        data_input <- read.csv(file_name, check.names = "false")
    } else if (file_extension %in% c("tsv", "txt", "tabular")) {
        data_input <- read.delim(file_name, sep = "\t", check.names = "false")
    } else if (file_extension == "parquet") {
        data_input <- arrow::read_parquet(file_name)
    } else {
        stop("Unsupported file format.")
    }
    return(data_input)
}
