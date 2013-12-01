Galaxy wrapper for GATK2
========================

This wrapper is copyright 2013 by Björn Grüning, Jim Johnson & the Galaxy Team.

The Genome Analysis Toolkit or GATK is a software package developed at the 
Broad Institute to analyse next-generation resequencing data. The toolkit offers
a wide variety of tools, with a primary focus on variant discovery and 
genotyping as well as strong emphasis on data quality assurance. Its robust 
architecture, powerful processing engine and high-performance computing features
make it capable of taking on projects of any size.

http://www.broadinstitute.org/gatk
http://www.broadinstitute.org/gatk/about/citing-gatk


GATK is Free for academics, and fee for commercial use. Please study the GATK licensing website:
http://www.broadinstitute.org/gatk/about/#licensing


Installation
============

The recommended installation is by means of the toolshed_.

.. _toolshed: http://toolshed.g2.bx.psu.edu/view/bjoern-gruening/augustus

Galaxy should be able to automatically install samtools dependencies automatically
for you. GATK2, and its new licence model, does not allow us to distribute the GATK binaries.
As a consequence you need to install GATK2 by your own, please see the GATK website for more informations:

http://www.broadinstitute.org/gatk/download

Once you have installed GATK2 you need to edit the env.sh file that is installed with these wrappers.
You will find this env.sh file under:

<tool_dependency_dir>/gatk2/<version>/iuc/<hash_string>/env.sh

You should edit the GATK2_PATH environment variable to point to the folder you have installed GATK2
and if you want to deactivate the 'call home feature' from GATK you can set

GATK2_SITE_OPTIONS='-et "NO_ET" -K "/data/gatk2_key_file"'

GATK2_SITE_OPTIONS can be used to insert specific options into every GATK2 wrapper 
during runtime, without changing the actuall wrapper.

Read more about the "Phone Home" problem under:
http://www.broadinstitute.org/gatk/guide/article?id=1250


Finally, you should fill in additional information about your genomes and 
annotations in the gatk2_picard_index.loc and gatk2_annotations.txt. 
You can find them under ./tool-data/.



History
=======

v0.1 - Initial public release


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


