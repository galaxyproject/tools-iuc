EBI Search
==========

EBI Search is a tool to provide text search functionality and uniform access to resources and services hosted at the European Bioinformatics Institute.

As the possible options in EBI Search are numerous, the `macros.xml` for this wrapper with all options is automatically generated using [`ebeye_urllib3.py`](http://www.ebi.ac.uk/Tools/webservices/download_clients/python/urllib/ebeye_urllib3.py) tool from EBI and a Python script ([`generate_macros.py`](generate_macros.py)). 

For any change in the `macros.xml`, please change on [`generate_macros.py`](generate_macros.py) and regenerate the `macros.xml` with

```
$ conda env create -f environment.yml
$ source activate ebeye_urllib
(ebeye_urllib) $ python generate_macros.py
```