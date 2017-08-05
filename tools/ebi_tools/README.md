EBI Search
==========

EBISearch is a tool to provide text search functionality and uniform access to resources and services hosted at the European Bioinformatics Institute.

For any change in the `macros.xml`, please change on [`generate_macros.py`](generate_macros.py) and regenerate the `macros.xml` with

```
$ conda env create -f environment.yml
$ source activate ebeye_urllib
(ebeye_urllib) $ python generate_macros.py
```