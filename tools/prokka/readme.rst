Prokka wrapper
==============

Warning
-------

Prokka includes custom databases and is thus about a 340 MB download!

Version history
---------------

- Release 9 (prokka 1.12.0): Support Prokka 1.12. Uses updated conda dependency 1.12
- Release 8 (prokka 1.11.1): Use Conda for installing tool dependencies, ToolShed dependencies are not supported anymore. Add a test. The tool is now maintained by the IUC.
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
