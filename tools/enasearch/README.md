ENASearch
=========

[ENASearch](https://github.com/bebatut/enasearch) is a Python library for interacting with [ENA](http://www.ebi.ac.uk/ena/browse/programmatic-access)'s API.

For any change in the `macros.xml`, please change on [`generate_macros.py`](generate_macros.py) and regenerate the `macros.xml` with

```
$ conda install enasearch
$ python generate_macros.py
```