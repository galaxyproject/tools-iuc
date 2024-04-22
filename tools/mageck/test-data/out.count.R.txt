Sweave("output_countsummary.Rnw");
library(tools);

texi2dvi("output_countsummary.tex",pdf=TRUE);

