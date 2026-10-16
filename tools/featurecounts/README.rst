FeatureCounts wrapper for Galaxy
================================

* http://bioinf.wehi.edu.au/featureCounts/
* http://subread.sourceforge.net/

FeatureCounts as part of the SUBREAD package is "a highly efficient and
accurate read summarization program".

Installation
------------

This wrapper requires Galaxy 16.04 to be fully functional because
of the following commits:

* https://github.com/galaxyproject/galaxy/pull/961
* https://github.com/galaxyproject/galaxy/pull/1714

Built-in annotations
--------------------

The wrapper defines the four annotation keys shipped with Subread (``hg19``,
``hg38``, ``mm10``, ``mm9``) and filters them by the alignment's genome build.
No ``featurecounts_anno`` data table or external reference-data installation is
needed for built-in annotations. Existing installations can retain that table
for older wrapper versions; this version does not read it. The legacy registry
files and table configuration are retained for compatibility.

Reference-dependent tests
-------------------------

The built-in annotation tests explicitly select ``hg19`` and need no reference
data tables. Additional tests reject unsupported and mismatched genome builds.
The cached annotation test selects
``featurecounts_test_hg38`` from ``gene_sets``; this key identifies the small
bundled ``featureCounts_guide.gff`` fixture, not a production hg38 annotation.
``tool_data_table_conf.xml.test`` supplies the cached fixture table for isolated
tests and retains the legacy built-in registry configuration.

The explicit cached selection exposes its ``gene_sets`` dependency to clients
running tests against an existing Galaxy server. Such clients should apply
their reference-data policy or provide an isolated fixture table rather than
relying on the first production annotation with a matching genome build.

License
-------

**featureCounts**:

GPL (>=3)

