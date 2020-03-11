#!/usr/bin/env bash

# install conda
if type conda > /dev/null; then  
	true
else
	tmp=$(mktemp -d)
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p "$tmp/miniconda"
	source "$tmp/miniconda/bin/activate"
fi

eval "$(conda shell.bash hook)"

# install conda env
if grep -Fq __bioconductor-dada2@1.14 <<< $(conda env list); then
	true
else
	conda create -y --quiet --override-channels --channel conda-forge --channel bioconda --channel defaults --name __bioconductor-dada2@1.14 bioconductor-dada2=1.14
fi

conda activate __bioconductor-dada2@1.14

# create test data
cd test-data/

# download Mothur SOP data from zenodo (GTN), same as 
# http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip but stable links
# but file names need to be fixed
wget -nc -O F3D0_S188_L001_R1_001.fastq https://zenodo.org/record/800651/files/F3D0_R1.fastq?download=1
wget -nc -O F3D0_S188_L001_R2_001.fastq https://zenodo.org/record/800651/files/F3D0_R2.fastq?download=1
wget -nc -O F3D141_S207_L001_R1_001.fastq https://zenodo.org/record/800651/files/F3D141_R1.fastq?download=1
wget -nc -O F3D141_S207_L001_R2_001.fastq https://zenodo.org/record/800651/files/F3D141_R2.fastq?download=1

# zip and reduce data to ~ 10% (for speed)
for i in *fastq
do
	head -n 3000 "$i" | gzip -c > "$i.gz"
done
rm *fastq

# download data bases from https://zenodo.org/record/158955
# as mentioned in https://benjjneb.github.io/dada2/training.html
wget -nc -O reference.fa.gz https://zenodo.org/record/158955/files/rdp_train_set_14.fa.gz?download=1
wget -nc -O reference_species.fa.gz https://zenodo.org/record/158955/files/rdp_species_assignment_14.fa.gz?download=1

# take ~ 5% of the reference (for speed)
zcat reference.fa.gz | head -n 1000 | gzip -c > t && mv t reference.fa.gz
zcat reference_species.fa.gz | head -n 1000 | gzip -c > t && mv t reference_species.fa.gz


# generate outputs
Rscript gentest.R

conda deactivate

# # remove files only needed for test generation
# rm learnErrors_F3D0_R2.pdf dada_F3D0_R2.Rdata
