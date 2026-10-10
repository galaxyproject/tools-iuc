Savont
======

Wrappers for `Savont <https://github.com/bluenote-1577/savont>`_, a tool that
turns high-accuracy long-read amplicon sequences (ONT R10.4.1 or PacBio HiFi)
into amplicon sequence variants (ASVs) at single-nucleotide resolution and
profiles their taxonomy against reference databases.

- Savont ASV (``savont asv``)
- Savont classify (``savont classify``)
- Savont SINTAX (``savont sintax``)
- Savont export (``savont export``)

Savont is available on `Bioconda <https://bioconda.github.io/recipes/savont/README.html>`_.

Important
---------

Savont expects full-length amplicon reads with a base accuracy of at least
98% (e.g. ONT R10.4.1 or PacBio HiFi). Standard short-read (Illumina) data is
**not** a suitable input.

Reference databases
-------------------

The **Savont classify** and **Savont SINTAX** tools classify against Savont
reference databases (EMU, SILVA, GreenGenes2, UNITE). Databases are installed
with the *Savont reference databases* data manager
(`data_manager_savont_dbs <https://github.com/galaxyproject/tools-iuc/tree/main/data_managers/data_manager_savont_dbs>`_).

Typical workflow
----------------

1. **Savont ASV**: generate ASVs and a feature table from the raw reads.
2. **Savont classify** (species-level, alignment-based) or **Savont SINTAX**
   (genus-level, k-mer bootstrapping): assign a taxonomy to the ASVs.
3. **Savont export** *(optional)*: merge the results of several samples into
   a QIIME2-compatible feature table, representative sequences, taxonomy and
   taxon counts.

Test data
---------

The test data consists of small subsets (500 and 150 reads) of the ZymoBIOMICS
Microbial Community II standard reads provided with the Savont test suite, and
a mock EMU-format reference database built from the Savont reference ASVs, so
that tests can verify the database handling without downloading a multi-gigabyte
reference database.
