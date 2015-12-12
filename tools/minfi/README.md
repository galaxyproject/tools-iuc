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

The minfi analysis pipeline for TCGA data provides only two functions as of now, to find differentially methylated regions and differentially methylated positions. This is because the data which is provided/used for this tool has been put through quality control measures. The Level 3 data is ready to use for analysis. 

* The first step of this process will be to fetch TCGA data from http://gdac.broadinstitute.org/runs/ , where we would need Standard data from a particular date, or the latest which is available at http://gdac.broadinstitute.org/runs/stddata__latest/ . We want the "Open" data avaiable to download.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.42.17%20PM.png" width="200">

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.42.42%20PM.png">


* Choose the cancer type with the Illumina 450K methylation data, for the sake of this example we will use "UCEC" - Uterine Corpus Endometrial Carcinoma. In the directory structure, we want the data with the suffix ```__Level_3__within_bioassay_data_set_function__data.Level_3```. This file is usually the largest in the directory and can be easily spotted by sorting the directory based on Size. So for the example, we will choose the file ```gdac.broadinstitute.org_UCEC.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0.tar.gz```. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.41.05%20PM.png">

* There are two options to obtain the data, 1. Get the link to that file, and paste it in your Galaxy Upload tool to fetch the data. 2. Download the data on to your local machine, and upload the whole file.

Fetch the data at the URL directly to galaxy.
<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.30.15%20PM.png">

Upload your local file: 
<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.34.30%20PM.png">

* Once the upload is finished, it should be available in the History panel of the session. Choose the tool "Minfi Analysis pipeline for TCGA data". Galaxy at one step uncompresses your dataset to show you only a ".tar" file, choose that file as the input to your tool. 



* This example will show the tool run, with the Basic Default settings. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.34.53%20PM.png">

* After the tool is executed, the tool will show two results in the history panel.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.59.32%20PM.png" width="150">

* This shows the two results "Differentially methylated Regions", and "Differentially methylated Positions". You can ```View``` the results by clicking on the ```eye``` icon. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.59.59%20PM.png">


##Test Data

The test_data folder contains some test files which are,

1. 572364052 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

2. 572364053 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

3. RGset.RData : RGChannelSet in RData format, containing Illumina slides 572364052 and 572364053,
divided by phenotype (Case and Control).

4. example.tar.gz : Sample archived and compressed file containing two dummy samples. This is to test if the galaxy upload feature unarchives files.

5. bumps.csv : Differentially methylated regions as output.
