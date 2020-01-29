
conda create -y --quiet --override-channels --channel conda-forge --channel bioconda --channel defaults --name __bioconductor-dada2@1.12 bioconductor-dada2=1.12
conda activate bioconductor-dada2@1.12

Rscript gentest.R

# remove files only needed for test generation
rm learnErrors_F3D0_R2.pdf dada_F3D0_R2.Rdata
