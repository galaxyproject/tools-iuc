# Multi-align


## Description


Align independently sequences of a multi-fasta file

## Prerequisites


1. Unix system with A Galaxy server
2. Some tools and module are used by Multi-align and must be installed via conda dependancies, conda manage dependencies by search package specify by requirement tag in multi-align.xml

requirement for perl-data-dumper:

- perl-test-more v1.001002

requirement for perl-extutils-makemaker :
- perl-app-cpanminus v1.7043
- perl-data-dumper v2.161

requirement for perl-sys-cpu :

- perl-extutils-makemaker v7.24
- perl-threaded v5.22.0


tool requirement : 

- perl-getopt-long v2.50 
- perl-statistics-r v0.34
- perl-sys-cpu v0.61
- perl-parallel-forkmanager v1.17
- perl-string-random v0.30
- samtools v1.8
- bwa v0.7.17
- r-ggplot2 v2.2.1
- bioconductor-gviz v1.22.3
- python v3.6


all dependencies must be in perl 5.22

## Manual Installation


The process has to be completed by an administrator of your Galaxy server to install  Multi-align

1. make a folder with multi-align.pl,multi-align.xml,RNAseq.pm,align_count_stranded.pl,G_windows.pl and bootstrap folder.

2. Put the tool into Galaxy's tools directory You need to add files into tools/ directory , where all tool-related files are stored, within your Galaxy installation.
3. 
Make Galaxy aware of the new tool Now that the tool and its definition file are ready, the final step is to make Galaxy aware of the new files. Galaxy recognizes installed tools by reading the tool_conf.xml tool configuration file. Thus, letting Galaxy know about the new tool is as easy as adding a few lines to the tool_conf.xml file located in the config/ directory of the Galaxy installation. New tools can either be added to existing sections or added to new sections defined in the following way:
    
    ```xml
     <section name="NewTools" id="mTools">
        <tool file="multi-align.xml" />
     </section>
    ```
  
4. Start or Restart Galaxy to use it. 
5. Add tool to the whitelist (from admin panel) to bypass sanitization.
