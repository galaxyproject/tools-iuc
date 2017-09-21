========================
Galaxy multiqc wrapper
========================

Aggregate results from bioinformatics analyses across many samples into a single report

MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.

============
Installation
============

Requirements: pkg-config

Conda recipe: A conda recipe exists in bioconda repository. If you have a recent Galaxy release, the packages installation is not mandatory.

================
Prevent displaying MultiQC webpage as gibberish
================

For Galaxy to display MultiQC's HTML output properly, you need to either:

1. Deactivate the sanitize_all_html option in galaxy.ini (sanitize_all_html = False), or
2. Whitelist the tool in "Manage Display Whitelist" after installing

=======
History
=======

 * v1.0.0.0	Update to multiQC 1.0, add additional supported tools
 * v0.6:        Initial public release

==========
Citation
==========

Created by Enancio (originally EnginesOn start_up project) company
