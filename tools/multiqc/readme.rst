Galaxy MultiQC wrapper
========================

Aggregate results from bioinformatics analyses across many samples into a single report

MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.


Prevent displaying MultiQC webpage as gibberish
-----------------------------------------------

For Galaxy to display MultiQC's HTML output properly, you need to either:

1. Deactivate the sanitize_all_html option in galaxy.ini (sanitize_all_html = False), or
2. Whitelist the tool in "Manage Display Whitelist" after installing

Support new modules
-------------------

Currently, the wrapper supports the modules for tools found on the MTS.
To add new ones, you can look at the patterns at https://github.com/ewels/MultiQC/blob/master/multiqc/utils/search_patterns.yaml
