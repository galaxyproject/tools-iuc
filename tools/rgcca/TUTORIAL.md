# TUTORIAL RGCCA GALAXY-TOOL 

##### Version: 1.0

##### Author: Etienne CAMENEN

##### Key-words: 
omics, RGCCA, multi-block

##### EDAM operation: 
analysis, correlation, visualisation

##### Contact: 
arthur.tenenhaus@l2s.centralesupelec.fr

##### Short description:
Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.) and produces textual and graphical outputs (e.g. variables and individuals plots).

---

## Description

We consider J data matrices X1 ,..., XJ. Each n × pj data matrix Xj = [ xj1, ..., xjpj ] is called a block and represents a set of pj variables observed on n individuals. The number and the nature of the variables may differ from one block to another, but the individuals must be the same across blocks. We assume that all variables are centered. The objective of RGCCA is to find, for each block, a weighted composite of variables (called block component) yj = Xj . aj, j = 1 ,..., J (where aj is a column-vector with pj elements) summarizing the relevant information between and within the blocks. The block components are obtained such that (i) block components explain well their own block and/or (ii) block components that are assumed to be connected are highly correlated. In addition, RGCCA integrates a variable selection procedure, called SGCCA, allowing the identification of the most relevant features (see [here](https://github.com/rgcca-factory/RGCCA/tree/release/3.0.0#algorithm) for more information).


## Setting up (not for Windows)

### Conda installation

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh [yes, then ./miniconda3 and no]
source ./miniconda3/bin/activate
```

### Package installation

```
PACKAGE_NAME=RGCCA
git clone https://github.com/rgcca-factory/${PACKAGE_NAME}.git
cd ${PACKAGE_NAME}
conda create -n ${PACKAGE_NAME} -c conda-forge -c bioconda -c icm-iconics planemo
```

### Usage
```
source ../miniconda3/bin/activate ${PACKAGE_NAME} [if environment is not activated]
planemo serve [wait, then go to the link]
```


## 1. Load the inputs


In the tool-shed (left panel), select the « RGCCA » tool (Fig. 1). 

Download the pre-formatted files [here](https://github.com/rgcca-factory/RGCCA/tree/release/3.0.0/inst/extdata). This folder includes three blocks with the same individuals (corresponding to the countries here) but different types of variables (agriculture, industry and politic). In this dataset, according to Russett (1964), a high agriculture inequality and a low industrial development lead to unstable political regime. 

Download them in Galaxy (with the download button in green, **Fig. 1**). The accepted format is one (for PCA) or multiple CSV files containing a matrix with:
- quantitative values only, with decimals separated by '.' and missing values labelled as "NA"
- samples in rows, labelled in the 1rst column with the same sample names between blocks (some samples could be missing in some blocks)
- variables in columns, labelled in the 1rst line without duplications in variable names between blocks

![Fig 1](static/images/toolShed.png)

*Fig. 1 : Tool-shed of Galaxy with the "RGCCA" emplacement*

The structure of the dataset should be seen in the history panel (**Fig. 2**).

![Fig 2](static/images/history.png) 

*Fig. 2. History of Galaxy after downloading the three blocks from Russett data*

The block could be then used in the "RGCCA" tool (**Fig. 3**).  
Use ```agriculture.tsv``` as dataset. Click on "Insert New dataset" to make new panels appear and to add 
```industry.tsv``` and ```politic.tsv``` as new dataset.

![Fig 3](static/images/tool.png) 

*Fig. 3. Graphical interface in Galaxy of "RGCCA". By default, an only parameter is required : an input file to analyze. Another dataset could be added for a multi-bloc analyze (in red).*


## 2. Customize the parsing


All the files should be delimited by the same type of separator (tabulation, by default) and by considering the first row as header. These parameters could be customized  by clicking on the "eye" icon in parsing settings and selecting semicolon (**Fig 4**). 

Change the separator to another one (e.g.,```semicolon```) to observe an error notification after executing the tool: 
> "politic block file has only one column. Check the separator."

![Fig 4](static/images/advParse.png) 

*Fig. 4. The panel of parsing settings should appear after clicking on the "eye" icon*


## 3. Custom the analysis


The analyse parameters are all set by default (**Fig. 5**) and the user could directly click on the ```Execute``` button. To directly visualize the outputs, see the [last section](https://github.com/BrainAndSpineInstitute/rgcca_galaxy/tree/master#5-visualize-the-plot).

![Fig 5](static/images/advAn.png) 

*Fig. 5 : The second parameter panel shows various options to customize the analysis: choose the analysis and the number of components, scale the blocks, choose a shrinkage, use the superblock or a supervised approach, choose a link function.*

### 3.1. Number of components and scaling

With all analysis methods, the ```number of components``` could be changed. By default, it is set to two components (for biplots). In the software, the maximum of components allowed is *a posteriori* limited by the minimum number of columns between all datasets. In the case of Russet data, two components are allowed because of the two columns in the industry block. Select only the agriculture and the politic bloc to move the cursor to three components. In any cas, five components are the maximum value allowed.

One could also selects ```scale [/ unscale] the blocks```. Either the option is selected or not, a data centering step is always performed. If selected, each block is normalized and then divided by the square root of its number of variables. When the data are already scaled, this step could be avoided by disabling the parameter.

### 3.2. Analysis methods

By default, the selected ```analysis method``` is set on ```RGCCA```. This tutorial will be focused on the RGCCA case, but another methods could be selected. For example, when only one block file is loaded in the previous step, a ```PCA``` could be performed. 

### 3.3. Connection between blocks

This parameters are only used only with R/SGCCA.

#### 3.3.1. Loading a connection file

The downloaded folder contains a design matrix (```connection.tsv```; **Fig. 6**) corresponding to the relationship between each block: 1 if two blocks are connected and 0 otherwise. The expected format should be tabulation-separated and do not have column and row names. It is a symmetric matrix with the same dimension as the number of blocks. This file allows to add *a priori* information of correlation hypothesis between the blocks. It will not be taken in account with a superblock (see next section). After disabling the ```use a superblock``` option, load this file into the ```Connection design``` file box and observe the result on the plots. The ```connection.tsv``` file contains 1 in all non-diagonal cells and makes the assumption that all the blocks are related.

![Fig 6](static/images/connection.png)

*Fig. 6. Supplementary files that could be used to customize the connection between the blocks in the analysis.*

#### 3.3.2. Superblock 
By default, all the blocks are connected to a superblock, a concatenation of all the other blocks. The space spanned by global components is viewed as a compromise space that integrated all the modalities and facilitates the visualization of the results and their interpretation. To visualize the blocks without the superblock, disable the ```Use a superblock``` option.

#### 3.3.3. Supervised analysis
By selecting ```supervised analysis``` option, a slider bar appears to select the block used as a response. By selecting this block, all other blocks (predictors) will be only connected to this block. For example, select the ```1``` for the agriculture block (the first block loaded previously).

If a superblock is used, supervised analysis should be disabled, and inversely.

### 3.4. Other R/SGCCA parameters

#### 3.4.1. Shrinkage parameter (Tau)
```Optimal tau``` is automatically set for each block. When disabled, one could make ```tau``` varying for each block from 1 (maximize the correlation between the variables of the selected block) to 0 (maximize the covariance).

#### 3.4.2. Sparsity coefficient
By selecting a SGCCA, the ``` sparsity``` could be applied to each block. This coefficient varies from the inverse of the square root of the number of columns (the smaller set of variables) to 1 (all the variables are included).

After selecting the SGCCA mode, move  the cursor to the ```0.75``` value to select less variables.

#### 3.4.3. Scheme function (advanced users)
```Scheme function``` allows to select the link (i.e. scheme) function for covariance maximizations between block components among: 
- identity (```Horst```)
- absolute values (```centroid```)
- squared values (```factorial```)

Only, the horst scheme penalizes structural negative correlation. The factorial scheme discriminates more strongly the blocks than the centroid one.


## 4. Customize the graphics


![Fig 7](static/images/advGraph.png) 

*Fig. 7 : The graphical option panel includes: (i) the loading of groups of response to color the samples (ii) the possibility to hide/print the names of the variables, (iii) the components used in the plots and (iv) the selection of the block to visualize. In this example, the superblock will be selected as the block for the X- and Y-axis.*

### 4.1. Color the samples
A variable could be used to color the points according to a response. For this, load the ``` political_system.tsv``` file (**Fig. 8**) in the corresponding ```groups of modalities``` box to update the plot. The expected format is a CSV file tabulation-separated with: 
- qualitative or quantitative values (decimals separated by '.') with missing values labelled as "NA"
- samples in lines, labelled in the 1rst column with the same sample names as the blocks (some samples could be missing)
- a header containing the names of the columns

![Fig 8](static/images/response.png)

*Fig. 8. Supplementary files that could be used to visualize the group of a response in the samples plot.*

### 4.2. Display names
If activated (by default), the ```display names``` option shows the name of the points in the biplots. If disabled, shapes are shown instead of text: one per group of modality.

### 4.3. Components (for the x/y-axis)
The ```component``` of the analysis allows to choose the space where the points are visualised. For the "top variable" histogram, the component is set by ```component for the x-axis```. For biplots tabs, either ```component for the x-axis``` or ```component for the y-axis``` could be set. By default, they are respectively set to the first and the second components. Their choices are limited by the number of components selected in the analysis (defined in the 2.2. section). If the number of components in RGCCA were greater than two (not allowed in the Russet example, because of the industry block), the ```component for the x-axis```, for example, could be set to the third one.

### 4.4. Block (for the x/y-axis)
By default, plots are shown with the ```superblock``` (i.e., the concatenation of all blocs; see [section 3.3.2](https://github.com/BrainAndSpineInstitute/rgcca_galaxy/tree/master#332-superblock)) to visualize all the blocs together. If this option is disabled, by default, the last blocks in the drop-down menu ```block``` is used (option ```0```). Choose another block (e.g., ```1``` for agriculture) to update the plots with your selection. For the ```individual``` plot, a ```block for the y-axis``` could also be selected.

### 4.5. Number of top variables
Used only in the  the "top variable" histogram, the maximum ```number of top variables``` is *a posteriori* automatically set to the number of variables in the selected blocks. For example, with Russet data, eleven "top" variables could be visualised by default on the superblock.


## 5. Visualize the plot 


Please, make sure that these options are set to visualise the same plots than those in the next examples:
- RGCCA method
- two components
- block scaled
- optimal tau
- a superblock
- factorial scheme

By executing the analysis (blue button at the bottom), four images, two tabular files and a RData should appear in the history panel. For each axis of the block, the corresponding percent of average explained variance is indicated in corresponding images.


### 5.1. Connection between blocks
```design.pdf``` summarizes the connection between each block: a link corresponds to a "1" value, in the matrix connection file (**Fig. 9**; see [section 3.3.1.](https://github.com/BrainAndSpineInstitute/rgcca_galaxy/tree/master#331-loading-a-connection-file)). For each block:
- "P" is the number of variables
- "N" is the number of lines (here, each block has the same number of line)
- "tau" is the shrinkage parameter and "sparsity" is the sparsity coefficient (see the [3.4.1 & 3.4.2 sections](https://github.com/BrainAndSpineInstitute/rgcca_galaxy/tree/master#34-other-rsgcca-parameters). The tau parameter could be shown for each component if the optimal option is selected

![Fig 9](https://raw.githubusercontent.com/rgcca-factory/RGCCA/release/3.0.0/img/design.png)

*Fig. 9 : Connection between each block of the RGCCA and the superblock with 47 common rows between blocks*

### 5.2. Average variance explained (AVE)
In ```ave.pdf``` the average variance explained (AVE; in X-axis) is represented in percent for each block (in Y-axis) and each component (one color per component) (**Fig. 10**). The subtitle informs about the AVE for the two first of the outer model (weighted average of the AVE of each block).

![Fig 10](https://raw.githubusercontent.com/rgcca-factory/RGCCA/release/3.0.0/img/ave.png)

*Fig. 10 : Average variance variance explained (in %) for each block and for the two first components of the RGCCA*

### 5.3. Samples
```individual.pdf``` is the projection of the sample coordinates in the selected component of the analysis and, by default, on the superblock (a concatenation of all the blocks) (**Fig. 11**). If a ```response``` file is loaded, each sample is colored according to this variable. In the Russet example, the X-axis could discriminate a dictatorship (with upper values on this axis than the two other political systems), whereas the Y axis discriminates an unstable democracy (with upper values than the others).

![Fig 11](https://raw.githubusercontent.com/rgcca-factory/RGCCA/release/3.0.0/img/individuals.png)

*Fig. 11 : Samples coordinates on the two first components for the superblock of the RGCCA by loading the " political_system.tsv" file.*

### 5.4. Corcircle
```corcircle.pdf``` corresponds to the Pearson correlation between the variables of the block and the selected components in the analysis (by default, on the two first components) (**Fig. 12**). The circle is a 1 correlation and the dotted one is a 0.5 correlation. If the superblock is selected, colors correspond to the belonging of each variable to each block. Only the 100th variables the most correlated to each axis are printed.

![Fig 12](https://raw.githubusercontent.com/rgcca-factory/RGCCA/release/3.0.0/img/corcircle.png)

*Fig. 12 : Correlation between each variable of each block (by using the superblock) and the two first components of RGCCA*

### 5.5. Top variables
```top_variables``` also represents the same correlation of the variable with the selected component (on the X-axis; 1 by default). The top variables are sorted decreasingly (on the Y-axis) in a histogram among the selected block (superblock, by default) (**Fig. 13**). Their number is set by the graphical associated parameter.

![Fig 13](https://raw.githubusercontent.com/rgcca-factory/RGCCA/release/3.0.0/img/top_variables.png)

*Fig. 13 : Top 11 variables among all the blocks (by using the superblock) with higher correlation with the first component of the RGCCA. "Gnpr" (belonging to the industry block) shows a correlation of 0.859 with this component.*

Here, "labo" from the industry block (the amount of labor force in agriculture) is the variable the most positively correlated to the X-axis. On the opposite, "gnpr" (gross national product) from the industry block and "demostab" (stable democracy) from the political block are the most negatively correlated variables. In other terms, these variables are the most importants on the first component of the RGCCA. Countries with an unstable democracy are more associated with a lower "rent" (percent of farmers that rent their land). Otherwise, those with a dictatorship system are more associated with a higher labor force in agriculture values and a less gross national product (and inversely for the stable democracy case).