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

Reference-dependent tests
-------------------------

The built-in annotation tests explicitly select ``hg19`` from
``featurecounts_anno``. The cached annotation test selects
``featurecounts_test_hg38`` from ``gene_sets``; this key identifies the small
bundled ``featureCounts_guide.gff`` fixture, not a production hg38 annotation.
``tool_data_table_conf.xml.test`` supplies these tables for isolated tests.

Explicit selections also expose the required data tables to clients running
the tests against an existing Galaxy server. Such clients should apply their
reference-data policy or provide isolated fixture tables rather than relying
on the first production annotation with a matching genome build.

License
-------

**featureCounts**:

GPL (>=3)

