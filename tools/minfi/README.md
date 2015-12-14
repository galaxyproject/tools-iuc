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

https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html

##Functions provided and description

1. **Minfi** **Pipeline** to analyze illumina 450k data. (minfi_pipeline.xml)

2. **Minfi** **Analysis** **pipeline** for TCGA data hosted on GDAC-Broad Institute.(minfi_TCGA_pipeline.xml) 


##Minfi Pipeline to analyze illumina 450k data 

####Tool xml name: *minfi_pipeline.xml* 

####How to use:


* The first step is to upload Illumina 450k IDAT files using the Galaxy Upload feature.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.43.07%20PM.png" >

* In this example, we uploaded "Case" samples from slide 572364052 - one Illumina slide with 6 IDAT files (3 samples with red/green channel pairs) - and "Control" samples from slide 572364053 - one Illumina slide with 6 IDAT files (3 samples with red/green channel pairs) - for a total of 12 IDAT files. Once uploaded, these files should appear in the history panel. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.43.41%20PM.png">

* The next step is to group samples based on treatment type. First, select "Operations on multiple datasets" in your history panel.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.44.32%20PM.png" width="150">

* Next, select the samples which belong to the "Case" treatment group (572364052_*.idat), and "For all selected", "Build a list of dataset pairs". 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.45.18%20PM.png" width="150">

* The "Create a collection of paired datasets" window will appear with the samples you selected. To pair the red/green channels, first "Clear Filters" and then type each dataset suffix - "_Grn.idat" and "_Red.idat" - in the filter bars. This should filter your list and display the "Green" channel files on one side and "Red" channel files on the other. If the datasets are correctly paied, click "Pair these datasets" in the middle of the window for each pair. Finally, name the collection - in this example "Case" based on the treatement type - and click "Create list".

![pairCollection](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.46.59%20PM.png)

* Repeat these steps to build a paired dataset for the "Control" samples (572364053_*.idat).

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.47.30%20PM.png" width="150">

* In addition to the option to pair each dataset, the "Auto-pair" link will pair all datasets at once. Again, after pairing, name the collection - in this example "Control" - and click "Create list".

![pairControl2](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.47.59%20PM.png)

* The history panel now displays two paired dataset collections: 13. "Case" and 14. "Control". 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.48.22%20PM.png" width="150">

* From the "Minfi" tools, choose "Minfi pipeline to analyse Illumina 450k data". In the tool form, select "Case" for "Condition1/Treatment" and "Control" for "Condition2/Wildtype". After choosing the "Quantile Normalization" preprocessing method and "Basic Default settings" for Minfi Parameters, run the tool.

![minfitoolform](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-08%20at%2010.30.54%20AM.png)

* Successful completion of the tool results in 4 items in the history panel: QC report, MDS plot, Differentially methylated positions, and Differentially methylated regions.

![result](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.53.03%20PM.png)

* To view results, click the "eye" in the history panel corresponding to the result.

![resultview](https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-07%20at%205.55.57%20PM.png)

##Minfi Analysis Pipeline for TCGA data hosted on GDAC Broad Institute 

####Tool xml name: *minfi_TCGA_pipeline.xml* 

####How to use:

The minfi analysis pipeline for TCGA data provides only two functions as of now, to find differentially methylated regions and differentially methylated positions. This is because the data which is provided/used for this tool has been put through quality control measures. The Level 3 data is ready to use for analysis. 

* The first step of this process is to fetch TCGA data from http://gdac.broadinstitute.org/runs/ where we need Standard data from a particular date, or the latest data which is available at http://gdac.broadinstitute.org/runs/stddata__latest/. We want the "Open" data avaiable to download.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.42.17%20PM.png" width="300">

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.42.42%20PM.png">


* Choose the cancer type with the Illumina 450K methylation data, for the sake of this example we will use "UCEC" - Uterine Corpus Endometrial Carcinoma. In the directory structure, we want the data with the suffix ```__Level_3__within_bioassay_data_set_function__data.Level_3```. This file is usually the largest in the directory and can be easily spotted by sorting the directory based on Size. So for the example, we will choose the file ```gdac.broadinstitute.org_UCEC.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0.tar.gz```. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-09%20at%202.41.05%20PM.png">

* There are two options to obtain the data: 1. Get the link to that file and paste it in your Galaxy Upload tool to fetch the data; 2. Download the data on to your local machine and upload the whole file.

Option 1 : Fetch the data at the URL directly to Galaxy

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.30.15%20PM.png">

Option 2: Upload a local file

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.34.30%20PM.png">

* Once the upload is finished, it should be available in the History panel of the session. Choose the tool "Minfi Analysis pipeline for TCGA data". 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-12%20at%205.17.44%20PM.png">


* Galaxy will uncompress your dataset to show you only a ".tar" file. Choose that file as the input to the tool. This example will show the tool run with the Basic Default settings.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.34.53%20PM.png">

* After the tool is executed, two results will appear in the history panel.

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.59.32%20PM.png" width="150">

* This shows the two results "Differentially methylated Regions" and "Differentially methylated Positions". You can ```View``` the results by clicking on the ```eye``` icon. 

<img src="https://github.com/nitesh1989/tools-iuc/blob/methylation_2/tools/minfi/help/help-images/Screen%20Shot%202015-12-11%20at%204.59.59%20PM.png">


##Test Data

The test_data folder contains some test files which are,

1. 572364052 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

2. 572364053 : One Illumina slide with 6 IDAT files (3 samples with red/green channel pairs)

3. dmrs.csv : Differentially methylated regions as output.
