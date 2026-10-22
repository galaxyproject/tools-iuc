# VirHunter

VirHunter is a deep learning method that uses Convolutional Neural Networks (CNNs) and a Random Forest Classifier to identify viruses in sequening datasets. More precisely, VirHunter classifies previously assembled contigs as viral, host and bacterial (contamination).

## System Requirements
VirHunter installation requires a Unix environment with [python 3.8](http://www.python.org/). 
It was tested on Linux and macOS operating systems. 
For now, VirHunter is still not fully compatible with M1 chip MacBook.

In order to run VirHunter you need to have git and conda already installed. 
If you are installing conda for the first time, we suggest you to use 
a lightweight [miniconda](https://docs.conda.io/en/latest/miniconda.html).
Otherwise, you can use pip for the dependencies' installation.
         
## Installation 

To install VirHunter, you need to download it from github and then to install the dependancies.

First, clone the repository from [github](https://github.com/cbib/virhunter)

git clone https://github.com/cbib/virhunter.git

Go to the VirHunter root folder

cd virhunter/

### Installing dependencies with Conda

First, you have to create the environment from the envs/environment.yml file.
The installation may take around 500 Mb of drive space. 

conda env create -f envs/environment.yml

Second, activate the environment:

conda activate virhunter

### Installing dependencies with pip

If you don't have Conda installed in your system, you can install python dependencies via pip program:

pip install -r envs/requirements.txt

Then if you have macOS you will need to install wget library to run some scripts (Conda installation already has it). You can do this with brew package manager.

brew install wget

### Testing your installation of VirHunter

You can test that VirHunter was successfully installed on the toy dataset we provide. 
IMPORTANT: the toy dataset is intended only to test that VirHunter has been well installed and all the scripts can be executed. 
These modules should not be used for prediction on your owd datasets!

First, you have to download the toy dataset

bash scripts/download_test_installation.sh

Then run the bash script that calls the testing, training and prediction python scripts of VirHunter.
Attention, the training process may take some time (up to an hour).

bash scripts/test_installation.sh


## Using VirHunter for prediction

To run VirHunter you can use the already pre-trained models or train VirHunter yourself (described in the next section).
Pre-trained model weights are already available for the multiple host plants. 
You can download them using the download_weights.sh script.

bash scripts/download_weights.sh

Once the config file is ready, you can start the prediction:

python virhunter/predict.py --test_ds /path/to/test_ds_1

After prediction VirHunter produces two csv files and one optional fasta file:

1. The first file ends with _predicted_fragments.csv
It is an intermediate result containing predictions of the three CNN networks (probabilities of belonging to each of the virus/plant/bacteria class) and of the RF classifier for each fragment of every contig.

2. The second file ends with _predicted.csv.
This file contains final predictions for contigs calculated from the previous file. 
   - id - fasta header of a contig.
   - length - length of the contig.
   - # viral fragments, # plant fragments and # bacterial fragments - the number of fragments of the contig that received corresponding class prediction by the RF classifier.
   - decision - class given by the VirHunter to the contig.
   - # viral / # total - number of viral fragments divided by the total number of fragments of the contig.
   - # viral / # total * length - number of viral fragments divided by the total number of fragments of the contig multiplied by contig length. It is used to display the most relevant contigs first.

3. The fasta file ends with _viral.fasta. It contains contigs that were predicted as viral by VirHunter.
