[![Build Status](https://travis-ci.org/kpbioteam/ewas_galaxy.svg?branch=master)](https://travis-ci.org/kpbioteam/ewas_galaxy)


# 🔬📚 Infinium Human Methylation BeadChip: a pipeline for population epigenetics integrated into Galaxy

This repository contains Infinium Human Methylation BeadChip pipeline that can be installed and used inside the Galaxy. 
The Tool is avaliable through [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/view/kpbioteam/ewastools/53aaf097238c) and ready to use via [Docker](https://galaxyproject.org/use/ewas-galaxy/).

Analysis pipline includes:
 * raw intensity data loading
 * .idat preprocessing
 * optional normalisation of the data and quality control with additional sample check.
   
This gives the user the opportunity to perform any of these preparation and data cleaning steps, including a recommended 'genetic variation annotation step'. Finally, the generated  dataset can be used to hunt (find) differentially-methylated positions DMP and regions DMR with respect to a phenotype co-variate.
