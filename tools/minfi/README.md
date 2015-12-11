Minfi (for galaxy)
====

*Maintainer: IUC*

*Contact: nitesh1989 on github*

The minfi package provides tools for analyzing Illumina’s Methylation arrays, with a special
focus on the new 450k array for humans.

In this package we refer to differentially methylated positions (DMPs) by which we mean
a single genomic position that has a different methylation level in two different groups of
samples (or conditions). This is different from differentially methylated regions (DMRs)
which imply more that more than one methylation positions are different between conditions.

### Table of Contents
**[User Guides](#user-guides)**  
**[Functions provided and description](#functions-provided-and-description)**  
**[Minfi Pipeline to analyze illumina 450k data](#minfi-pipeline-to-analyze-illumina-450k-data)**  
**[Minfi Analysis Pipeline for TCGA data hosted on GDAC Broad Institute](#minfi-analysis-pipeline-for-tcga-data-hosted-on-gdac-broad-institute)**  
**[Test Data](#test-data)**  


##User Guides
Here are some excellent user guides about minfi, written by the contributors.

http://www.bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.pdf

https://bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf

https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#bumphunter-to-find-differentially-methylated-regions-dmrs

##Functions provided and description

1. **Minfi** **Pipeline** to analyze illumina 450k data. (minfi_pipeline.xml)

2. **Minfi** **Analysis** **pipeline** for TCGA data hosted on GDAC-Broad Institute.(minfi_TCGA_pipeline.xml) 


##Minfi Pipeline to analyze illumina 450k data 

####Tool xml name: *minfi_pipeline.xml* 

####How to use:


* The first step is to upload your Illumina 450k IDAT files using the Galaxy Upload feature.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.43.07%20PM.png" >

* In the example being shown here, we include samples from slides 572364052 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs) and 572364053 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs). There are 12 IDAT files in total. Once these files are uploaded you should be able to see them in the history panel. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.43.41%20PM.png">

* You will now need to select operations on multiple datasets in your history panel.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.44.32%20PM.png" width="150">

* Select the samples which belong to a treatment type and "for all selected" , build a list of dataset pairs. In this step we are doing is grouping samples, based on the treatment type.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.45.18%20PM.png" width="150">

* You should now see the window "Create a collection of Paired Datasets" open up, with the samples you selected. First, "Clear Filters" and then type the suffix of your datasets, which should be "_Grn.idat" and "_Red.idat". This should filter your list, and show you all the "Green" channel on one side and "red" on the other. Now click "Pair these datasets" in the middle of the window for each pair. This is demonstrated in a clear manner in the following few screenshots. Now you should give your collection a name based on the treatement type, here we name it "Case".

![pairCollection](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.46.59%20PM.png)

* In your history panel now, you should build a list of dataset pairs for your "Control" samples (or second treatment type).

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.47.30%20PM.png" width="150">

* This shows the filtering of the datasets, and the option to "Pair" your IDAT files. After pairing, please name it based on a second treatement type, we call it "Control" in the example.

![pairControl2](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.47.59%20PM.png)

* Your history panel should now show, two sets of paired dataset collections, 13. "Case" and 14. "Control" in the history panel shown as an example. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.48.22%20PM.png" width="150">

* From the "Minfi" tools, choose "Minfi pipeline to analyse Illumina 450k data". In the tool form, we select, "Condition1/Treatment" as Case, and "Condition2/Wildtype" as Contorl. After choosing a preprocessing method "Qunatile" and Basic parameters, we run the tool.

![minfitoolform](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-08%20at%2010.30.54%20AM.png)

* This shows the successful run of the tool, with Basic parameters, and shows 4 results in your history panel. QC report, MDS plot, Differentially methylated positions and differentially methylated regions.

![result](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.53.03%20PM.png)

* To view the results, click on the "eye" in the history panel corresponding to the result. You should be able to see the output. 

![resultview](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.55.57%20PM.png)

##Minfi Analysis Pipeline for TCGA data hosted on GDAC Broad Institute 

####Tool xml name: *minfi_TCGA_pipeline.xml* 

####How to use:



##Test Data

The test_data folder contains some test files which are,

1. 572364052 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

2. 572364053 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

3. RGset.RData : RGChannelSet in RData format, containing Illumina slides 572364052 and 572364053,
divided by phenotype (Case and Control).

4. example.tar.gz : Sample archived and compressed file containing two dummy samples. This is to test if the galaxy upload feature unarchives files.

5. bumps.csv : Differentially methylated regions as output.
