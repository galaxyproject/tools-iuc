# package_minfi_1.15.7

Minfi (for galaxy)
====

The minfi package provides tools for analyzing Illumina’s Methylation arrays, with a special
focus on the new 450k array for humans.

In this package we refer to differentially methylated positions (DMPs) by which we mean
a single genomic position that has a different methylation level in two different groups of
samples (or conditions). This is different from differentially methylated regions (DMRs)
which imply more that more than one methylation positions are different between conditions.

User Guides:

http://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf

https://bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

Functions provided and description
===

1. Analyze Custom Illumina 450k IDAT files (minfi_pipeline.xml)

2. Analyze TCGA 450k IDAT files (minfi_TCGA_pipeline.xml) 

Steps to use *minfi_pipeline.xml*
===




Test Data
====

The test_data folder contains some test files which are,

1. 572364052 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

2. 572364053 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

3. RGset.RData : RGChannelSet in RData format, containing Illumina slides 572364052 and 572364053,
divided by phenotype (Case and Control).

4. example.tar.gz : Sample archived and compressed file containing two dummy samples. This is to test if the galaxy upload feature unarchives files.

5. bumps.csv : Differentially methylated regions as output.
