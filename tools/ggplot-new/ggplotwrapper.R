library(tidyverse)
library(rjson)
library(argparse)

parser <- ArgumentParser(description="ggplot2 wrapper")
parser$add_argument('--data', type="string", help="Input json string for ggplot")
args <- parser$parse_args("--data")

input <- fromJSON(args$data)

##READ IN FILE
file <- c(input$file,  input$header, input$rownames)
data <- read.table(file=input$file, header=input$header, rowm.names=input$rownames,  sep=input$sep)

##Create base ggplot command
aescount <-  0
commandstring <- paste0("ggplot(data, x=", input$x, ", y=", input$y, "aes(")
for (aesthetic in names(input$aes)){
  if(aescount != 0){
    commandstring <- paste0(commandstring, ",")
  }
  commandstring <- paste0(commandstring, aesthetic, "=", input$aes[[aesthetic]])
  aescount <- aescount + 1
}

commandstring <- paste0(commanstring, ")")


## Should now have "ggplot(data, x=*, y=*, aes(*+))" 
## Now add geometries

## Break into plots and aes
for (plot in input$plots){
  geomstring <- paste0(" + ", plot$plot_type, "(")
  aescount <- 0
  for(aesthetic in name(plot$aesthetics)){
    if (aescount != 0){
      geomstring <- paste0(geomstring, ",")
    }
    geomstring <- paste0(geomstring, aesthetic, "=", plot$aesthetics[[aesthetic]])
    aescount <- aescount + 1
  }
  geomstring <- paste0(geomstring, ")")
}

commandstring <- paste0(commandstring, geomstring)
outplot <- eval(commandstring)
outfile <- paste0("graph.", input$outformat)

ggsave(outfile, outplot, device=input$outformat)
