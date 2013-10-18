Prokka wrapper
==============

Warning
-------

Prokka includes custom databases and is thus about a 2.0 GB download!

Dependencies of Prokka which needs to be installed separately
-------------------------------------------------------------

- Perl core modules: File\::Copy, FindBin, Getopt::Long, List::Util, Scalar::Util, Time::Piece, Time::Seconds;
- Perl modules: Bio::SeqIO from BioPerl_ >= 1.6.900, `XML::Simple`_;
- `GNU Parallel`_ >= 20130422 is required, but is shipped with Prokka and thus is not managed by the tool dependency system;
- tbl2asn_ >= 21.0 is required. This dependency is not managed here since versions are increasing very rapidly;
- SignalP_ >= 3.0 is an optional dependency to find signal peptides. For licensing reasons, it is not used in the tool wrapper.

.. _BioPerl: http://search.cpan.org/dist/BioPerl/
.. _XML::Simple: http://search.cpan.org/dist/XML-Simple/
.. _GNU Parallel: http://www.gnu.org/software/parallel/
.. _tbl2asn: http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/
.. _SignalP: http://www.cbs.dtu.dk/services/SignalP/

Configuration
-------------

Change the PROKKA_SITE_OPTIONS variable in the installed env.sh file to adjust the number of CPUs to use (--cpus).

Version history
---------------

- v0.1 (LG): initial release in the toolshed, supports Prokka 1.6.
- v0.2 (LG): added this readme file, supports Prokka 1.7, and adds dependencies management.
- v1.1.0: merge the wrappers by CRS4 and Lionel Guy, add COPYING file, make all params optional, correctly quote text params in command, use float type for 'evalue' param, describe output files in help, upgrade BLAST+ dependency to version 2.2.28, depend on package_aragorn_1_2_36 instead of trna_prediction, add PROKKA_SITE_OPTIONS to env.sh and remove 'cpus' param.

