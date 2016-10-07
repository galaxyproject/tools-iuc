Prokka wrapper
==============

Warning
-------

Prokka includes custom databases and is thus about a 340 MB download!

Dependencies of Prokka which need to be installed separately
-------------------------------------------------------------

- Perl core modules: File\::Copy, FindBin, Getopt::Long, List::Util, Scalar::Util, Time::Piece, Time::Seconds;
- Perl modules: Bio::SeqIO from BioPerl_ >= 1.6.900, `XML::Simple`_;
- SignalP_ >= 3.0 is an optional dependency to find signal peptides. For licensing reasons, it is not used in the tool wrapper.

.. _BioPerl: http://search.cpan.org/dist/BioPerl/
.. _XML::Simple: http://search.cpan.org/dist/XML-Simple/
.. _SignalP: http://www.cbs.dtu.dk/services/SignalP/

Configuration
-------------

prokka tool may be configured to use more than one CPU core by selecting an appropriate destination for this tool in Galaxy job_conf.xml file (see https://wiki.galaxyproject.org/Admin/Config/Jobs and https://wiki.galaxyproject.org/Admin/Config/Performance/Cluster ).

If you are using Galaxy release_2013.11.04 or later, this tool will automatically use the number of CPU cores allocated by the job runner according to the configuration of the destination selected for this tool.

Version history
---------------

- Release 7 (prokka 1.11.0): Support Prokka 1.11. Upgrade dependencies to package_barrnap_0_7, package_blast_plus_2_2_31, package_hmmer_3_1b2, package_tbl2asn_24_3.
- Release 6 (prokka 1.4.0): Use <stdio> because prokka writes some warnings on stderr. Update Orione citation. Update Prokka citation. Support Prokka 1.10. Upgrade dependencies to package_minced_0_1_6, package_barrnap_0_5 and package_tbl2asn_23_7. Added --proteins option. Add <citations>.
- Release 5 (prokka 1.3.0): Fix Prokka 1.8 dependency installation.
- Release 4 (prokka 1.3.0): Support Prokka 1.8. Depend on package_minced_0_1_4 and package_tbl2asn_22_4 (requires Galaxy release_2013.11.04 or later). Update citation.
- Release 3 (prokka 1.2.0): Use $GALAXY_SLOTS instead of $PROKKA_SITE_OPTIONS. Upgrade Barrnap dependency to v. 0.3. Upgrade Infernal dependency to v. 1.1. Depend on package_gnu_parallel_20131122 (requires Galaxy release_2013.11.04 or later).
- Release 2 (prokka 1.1.0): Merge the wrappers by CRS4 and Lionel Guy. Directly call prokka, remove prokka.py . Add 'locustag', 'increment', 'gffver', 'compliant', 'addgenes', 'genus', 'species', 'strain', 'plasmid', 'gcode', 'usegenus', 'metagenome', 'fast', 'evalue', 'norrna', 'notrna' params. Upgrade BLAST+ dependency to v. 2.2.28. Add dependencies on prodigal and barrnap. Add readme.rst .
- Release 1 (prokka 1.0.1): Add txt output file. Use a definition list instead of a block quote in <help>. Correct 2 dependency minimum versions.
- Release 0 (prokka 1.0.0): Initial release in the Tool Shed.

Version history of (now deprecated) Lionel Guy's wrapper:

- prokka 1.1.0: Merge the wrappers by CRS4 and Lionel Guy. Add COPYING file with MIT license. Make all params optional. Add 'gffver' param. Correctly quote text params in command. Use float type for 'evalue' param. Describe output files in help. Upgrade BLAST+ dependency to v. 2.2.28. Depend on package_aragorn_1_2_36 instead of trna_prediction. Depend on package_prodigal_2_60 instead of prodigal. Depend on package_barrnap_0_2 instead of barrnap. Add PROKKA_SITE_OPTIONS to env.sh and remove 'cpus' param.
- prokka 0.2: Added this readme file. Support Prokka 1.7. Add dependencies management.
- prokka 0.1: Initial release in the Tool Shed, supports Prokka 1.6.

Development
-----------

Development is hosted at https://bitbucket.org/crs4/orione-tools . Contributions and bug reports are very welcome!
