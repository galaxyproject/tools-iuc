Wrappers for the core functionality of the dada2 package https://benjjneb.github.io/dada2/index.html. 

- filterAndTrim
- learnErrors
- dada
- mergePairs
- makeSequenceTable
- removeBimeraDenovo

Installation
============

The dada2 wrappers can be installed via the toolshed. Since they use datatypes that have been introduced with Galaxy release 19.09 they won't work out of the box for older Galaxy releases. 
In order to run the tools you may either upgrade Galaxy or execute the following two steps: 

1. `find GALAXY_ROOT/shed_tools/testtoolshed.g2.bx.psu.edu/repos/iuc/ -name "dada2_*xml" -exec sed -i -e 's/profile="19.09"/profile="YOUR_RELEASE"/' {} ;` (replace GALAXY_ROOT and YOUR_RELEASE appropriately)
2. insert the following lines in `config/datatypes.xml` (just before the line `</registration>`):
```
    <datatype extension="dada2_dada" type="galaxy.datatypes.binary:RData" subclass="true" display_in_upload="true" />
    <datatype extension="dada2_errorrates" type="galaxy.datatypes.binary:RData" subclass="true" display_in_upload="true" />
    <datatype extension="dada2_mergepairs" type="galaxy.datatypes.binary:RData" subclass="true" display_in_upload="true" />
    <datatype extension="dada2_sequencetable" type="galaxy.datatypes.tabular:Tabular" mimetype="application/text" subclass="true" display_in_upload="true" />
    <datatype extension="dada2_uniques" type="galaxy.datatypes.tabular:Tabular" mimetype="application/text" subclass="true" display_in_upload="true" />
```

Datatypes
=========

The dada2 Galaxy wrappers use a few extra data types to ensure that only inputs of the correct type can be used, these datatypes are available from Galaxy release 19.05, for earlier releases they need to be added manually. 

For the outputs of dada, learnErrors, and mergePairs the following datatypes are used that derive from  Rdata (which contains the named list that is returned from the corresponding dada function):

- dada2_dada (Rdata: named list, see docs for dada-class)
- dada2_errorrates (Rdata: named list, see docs for learnErrors)
- dada2_mergepairs (Rdata: named list, see docs for mergePairs)

For the outputs of makeSequenceTable and removeBimeraDenovo the following data types are used which derive from tabular:

- dada2_uniques
-- in R a named integer vector (names are the unique sequences)
-- in Galaxy written as a table (each row corresponding to a unique sequence, column 1: the sequence, column 2: the count)
- dada2_sequencetable
-- in R a named integer matrix (rows = samples, columns = unique sequences)
-- in Galaxy written as a table (rows = unique sequences, columns = samples)

Note the difference between the R and Galaxy representations! The main motivation is that the dada2_sequencetable is analogous to OTU tables as produced for instance by qiime (and it seemed natural to extend this to the uniques which are essentially a sequencetables of single samples).

Test data
=========

Test data for `dada2_seqCounts` is generated using planemo's `--update_test_data` argument and manual
inspection of the test files. In addition a run of the pipeline (using collections) is executed
manually using `planemo serve` making sure that the entries of the tables are generated in a useful way.

In order to have the Collection unzip tool available use `planemo s --galaxy_root GALAXY_ROOT  --extra_tools GALAXY_ROOT/lib/galaxy/tools/`

All test other test data is generated using the shell script (`gentest.sh`) in test-data 
