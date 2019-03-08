Berokka
=======

Trim, circularise, orient & filter long read bacterial genome assemblies

Introduction
------------

There is already a good piece of software to trim/circularise and orient
genome assemblies called `Circlator <https://sanger-pathogens.github.io/circlator/>`_.
Please try that first!

You should only try Berokka if:


1. You only have the contig files and do not have the corrected reads anymore
2. Your contigs are simple cases with clear overhang and could be done manually with BLAST
3. Circlator fails on your data even after `troubleshooting <https://github.com/sanger-pathogens/circlator/wiki/Troubleshooting>`_

**NOTE:** orientation to *dnaA* or *rep* genes is not yet implemented.

Installation
------------

Homebrew
^^^^^^^^

Using Homebrew will install all the dependencies for you:
`Linux <http://linuxbrew.sh>`_ or `MacOS <http://brew.sh>`_

.. code-block::

   brew install brewsci/bio/berokka

Conda
^^^^^

.. code-block::

   conda install -c bioconda  berokka

Source
^^^^^^

.. code-block::

   git clone https://github.com/tseemann/berokka.git
   ./berokka/bin/berokka -h

You will need to install all the dependencies manually:


* `BioPerl <http://bioperl.org/>`_ >= 1.6 (for ``Bio::SeqIO`` and ``Bio::SearchIO``\ )
* `BLAST+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ >= 2.3.0 (for ``blastn``\ )

Usage
-----

Input
^^^^^

Input should be completed long-read assemblies in FASTA format, such as those from
`CANU <https://github.com/marbl/canu>`_
or
`HGAP <https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP-in-SMRT-Analysis>`_.

Usage
^^^^^

.. code-block::

   % berokka --outdir trimdir canu.contigs.fasta
   <snip>
   Did you know? berokka is a play on the concept of overhang vs hangover

   % ls trimdir/
   01.input.fa
   02.trimmed.fa
   03.results.tab

   % cat trimdir/03.results.tab

   #sequence       status  old_len new_len trimmed
   tig00000000     trimmed 5461026 5448790 12236
   tig00000002     trimmed 138825  113601  25224
   tig00000003     trimmed 57075   43297   13778
   tig00000004     kept    24900   24900   0
   tig00000006     trimmed 1620    1320    300
   tig00000007     removed 2380    0       0

Output
^^^^^^

.. list-table::
   :header-rows: 1

   * - Filename
     - Format
     - Description
   * - 01.input.fa
     - FASTA
     - All the input sequences
   * - 02.trimmed.fa
     - FASTA
     - The (possibly) trimmed sequences
   * - 03.results.tab
     - TSV
     - Summary of results


The ``02.trimmed.fa`` output has been augmented with new header data (unless ``--noanno`` used):


* ``circular=true`` - inform that this is a circular sequence (Rebaler uses this)
* ``overhang=N`` - informs that N bp were trimmed off
* ``len=N`` - the new contig length if it was present (Canu adds this)
* ``suggestCircular=yes`` if the ``no`` version was present (Canu adds this)
* ``class=replicon`` if the ``class=contig`` was present *and* we circularised

Options
^^^^^^^


* 
  ``--filter <FASTA>`` allows you to remove contigs which match 50% of sequences in this file.
  Berokka comes with the standard Pacbio control sequence. You can provide your own FASTA file
  using this option. If you want to disable filtering, use ``--filter 0``.

* 
  ``--readlen LENGTH`` can be used for datasets that won't seem to circularise. 
  It affects the length of the match it attempts to make using BLAST.

* 
  ``--noanno`` will ensure that the FASTA descriptions are not altered between the input
  and output FASTA files.

* 
  ``--keepfiles`` and ``--debug`` are primarily for use by the developer.

Etymology
---------

`Berocca <https://en.wikipedia.org/wiki/Berocca>`_ is a brand of effervescent drink and vitamin tablets containing vitamin B and C.
It is a popular cure for a `hangover <https://en.wikipedia.org/wiki/Hangover>`_. A key role of the ``berokka`` tool is to remove the
"overhang" that occurs at the ends of long-read assemblies of circular genomes.

Feedback
--------

Please file questions, bugs or ideas to the `Issue Tracker <https://github.com/tseemann/berokka/issues>`_

License
-------

`GPLv3 <https://raw.githubusercontent.com/tseemann/berokka/master/LICENSE>`_

Citation
--------

Not published yet.

Authors
-------


* Torsten Seemann
