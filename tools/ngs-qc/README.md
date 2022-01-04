NGS-QC Generator
========

## About

Comparative analysis between ChIP-seq and other NGS datasets requires prior characterization of their degree of technical similarity.
NGS Q/C Generator infer Quality indicators from the distribution of sequenced reads associated to a particular NGS profile.
Such information is then used for comparative purposes and for defining strategies aiming to improve the quality of the sample-derived datasets.

## Authors

[Mendoza-Parra et. all](http://dx.doi.org/10.1093%2Fnar%2Fgkt829)

Improved and maintained by Matthias Blum and Pierre-Etienne Cholley.

## Requirements

* sort from the GNU coreutils
* Python 2.7 and the following packages:
    * numpy
    * scipy (requires `lapack-devel`, and `blas-devel`)
    * matplotlib
    * Pillow (requires `libzip-devel`, `libjpeg-devel`, and `libpng-devel`)
    * reportlab
* [Bedtools](http://bedtools.readthedocs.org/)
* [Gnuplot](http://www.gnuplot.info/) 4.4+ (requires `libcairo2-dev`, `libgd2-xpm-dev`, and `libpango1.0-dev`)

## Installation

Build the C tools by running `sh build.sh` in the **utils** directory.
Create a copy of `config.ini.sample`, rename it `config.ini`, and modify it (change the `softwares` section)


