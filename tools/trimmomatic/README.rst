Trimmomatic: flexible read trimming tool for Illumina NGS data
==============================================================

Galaxy tool wrapper for the Trimmomatic program, which provides various functions for
manipluating Illumina FASTQ files (both single and paired-end).

Trimmomatic has been developed within Bjorn Usadel's group at RWTH Aachen university
http://www.usadellab.org/cms/index.php?page=trimmomatic

The reference for Trimmomatic is:

- Bolger, A.M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer
  for Illumina Sequence Data. Bioinformatics, btu170.

Controlling the available memory
================================

The default amount of memory avilable to trimmomatic is set to 8GB.
To change the default amount of memory you can set the environment variable
``_JAVA_OPTIONS`` to ``-Xmx<amount_of_memory_in_GB>G``. The recommended way to
set this is in the job_conf.xml file. To change the available memory to 6GB, a
line like the below should be added:

``<env id="_JAVA_OPTIONS">-Xmx6G</env>``

This will set the environment variable ``_JAVA_OPTIONS`` to ``-Xmx6G``.

History
=======

============== ================================================================
Version        Changes
-------------- ----------------------------------------------------------------
0.39+galaxy1   - Relocated to the ``tools-iuc`` repository
0.39           - Update to Trimmomatic 0.39.
0.38.1         - Bug fix: add dependency on ``coreutils`` so that
                 ``readlink -e`` is supported across both Linux and MacOS
                 platforms.
0.38.0         - Update to Trimmomatic 0.38.
0.36.6         - Added trimlog and log outputs; add support for
                 ``fastqillumina`` and ``fastqsolexa`` input types
0.36.5         - Remove tool_dependencies.xml and always use conda to resolve
                 tool dependencies
0.36.4         - Add option to provide custom adapter sequences for
                 ILLUMINACLIP
               - Add options ``minAdapterLength`` and ``keepBothReads`` for
                 ILLUMINACLIP in palindrome mode
0.36.3         - Fix naming of output collections. Instead of all outputs being
                 called "Trimmomatic on collection NN" these will now be called
                 "Trimmomatic on collection NN: paired" or "Trimmomatic on
                 collection NN: unpaired".
0.36.2         - Support fastqsanger.gz datatype. If fastqsanger.gz is used as
                 input the output will also be fastqsanger.gz.
               - Use $_JAVA_OPTIONS to customize memory requirements.
0.36.1         - Reimplement to work with bioconda Trimmomatic 0.36 (toolshed
                 version is still supported for now).
0.36.0         - Update to Trimmomatic 0.36.
0.32.4         - Add support for ``AVGQUAL`` and ``MAXINFO`` operations.
0.32.3         - Add support for FASTQ R1/R2 pairs using dataset collections
                 (input can be dataset collection, in which case tool also
                 outputs dataset collections) and improve order and naming of
                 output files.
0.32.2         - Use ``GALAXY_SLOTS`` to set the appropriate number of threads
                 to use at runtime (default is 6).
0.32.1         - Remove ``trimmomatic_adapters.loc.sample`` and hard-code
                 adapter files into the XML wrapper.
0.32.0         - Add tool_dependencies.xml to install Trimmomatic 0.32
                 automatically and set the environment.
               - Update tool versioning to use Trimmomatic version number (i.e.
                 ``0.32``) with tool iteration appended (i.e. ``.1``).
0.0.4          - Specify '-threads 6' in <command> section.
0.0.3          - Added MINLEN, LEADING, TRAILING, CROP and HEADCROP options of
                 trimmomatic.
0.0.2          - Updated ILLUMINACLIP option to use standard adapter sequences
                 (requires the trimmomatic_adapters.loc file; sample version is
                 supplied) plus cosmetic updates to wording and help text for
                 some options.
0.0.1          - Initial version
============== ================================================================


Credits
=======

This wrapper was originally developed and maintained by Peter Briggs
(@pjbriggs).
Peter van Heusden (@pvanheus) and Marius van den Beek (@mvdbeek) contributed 
support for gz compressed FastQ files. Charles Girardot (@cgirardot) and
Jelle Scholtalbers (@scholtalbers) contributed additional options to
ILLUMINACLIP.
Matthias Bernt (@bernt-matthias) added log and trimlog output.
Nicola Soranzo (@nsoranzo) suggested using coreutils to enable cross-platform
support across Linux and MacOS.
Cristóbal Gallardo (@gallardoalba) updated Trimmomatic up to version 0.39.
Peter Briggs wishes to acknowledge the help from Matthia Bernt
(@bernt-matthias) with relocating the tool in the IUC tool repository,
and the IUC for taking on responsibility for the tool.

Developers
==========

The Trimmomatic tool is now maintained as part of the ``tools-iuc`` repository
on GitHub:
https://github.com/galaxyproject/tools-iuc/tools/tree/main/trimmomatic

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
