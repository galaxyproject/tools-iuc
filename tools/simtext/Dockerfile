FROM rocker/r-ubuntu:20.04
WORKDIR /simtext
COPY . .

RUN apt-get update -qq && apt-get -y upgrade && apt-get -y --no-install-recommends install \
  libssl-dev  \
  file \
  zlib1g-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  libjpeg-dev \
  && Rscript -e 'if (!require("shiny")) install.packages("shiny", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("plotly")) install.packages("plotly", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("DT")) install.packages("DT", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("shinycssloaders")) install.packages("shinycssloaders", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("shinythemes")) install.packages("shinythemes", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("tableHTML")) install.packages("tableHTML", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("argparse")) install.packages("argparse", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("PubMedWordcloud")) install.packages("PubMedWordcloud", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("ggplot2")) install.packages("ggplot2", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("stringr")) install.packages("stringr", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("tidyr")) install.packages("tidyr", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("magrittr")) install.packages("magrittr", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("plyr")) install.packages("plyr", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("ggpubr")) install.packages("ggpubr", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("rafalib")) install.packages("rafalib", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("RColorBrewer")) install.packages("RColorBrewer", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("dendextend")) install.packages("dendextend", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("Rtsne")) install.packages("Rtsne", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("umap")) install.packages("umap", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("reutils")) install.packages("reutils", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("textclean")) install.packages("textclean", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("easyPubMed")) install.packages("easyPubMed", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("RCurl")) install.packages("RCurl", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("SnowballC")) install.packages("SnowballC", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("mclust")) install.packages("mclust", repos="https://cloud.r-project.org");' \
  && Rscript -e 'if (!require("SemNetCleaner")) install.packages("SemNetCleaner", repos="https://cloud.r-project.org");' 

ENV PATH "/simtext/:${PATH}"
