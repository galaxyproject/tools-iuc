clinker
-------

This Galaxy wrapper runs `clinker`_ to compare annotated gene clusters or
plasmids and produce a standalone interactive HTML figure.

.. _clinker: https://github.com/gamcil/clinker

User-facing input, parameter, and output instructions are kept in the
``<help>`` section of ``clinker.xml``, where they are rendered in Galaxy.

Implementation notes
^^^^^^^^^^^^^^^^^^^^

The wrapper includes two helper scripts for behavior that is better handled in
Python than in the Galaxy command template:

* ``prepare_inputs.py`` validates every collection element as annotated
  GenBank, requires CDS features and at least two files, preserves collection
  order, and stages safe filenames for clinker. It also accepts the current
  Bakta ``annotation_gbff`` output when Galaxy labels it ``tabular``.
* ``sanitize_html.py`` parses clinker's embedded JSON and escapes
  HTML-sensitive characters before rewriting the output, preventing
  input-derived annotation text from closing the document's script element.

The scripts are listed in ``<required_files>`` so they are available when the
tool runs remotely.

Testing
^^^^^^^

Functional tests are defined in ``clinker.xml`` and run with Planemo. They
cover ordered collections, both GenBank and tabular-labelled GBFF inputs,
generic annotated records, invalid input, the one-file failure case, and a
regression test for HTML/JavaScript injection in annotation text.

License and citation
^^^^^^^^^^^^^^^^^^^^

The wrapper is distributed under the MIT license. Please cite:

Gilchrist, C. L. M. and Chooi, Y.-H. (2021). clinker & clustermap.js:
Automatic generation of gene cluster comparison figures.
https://doi.org/10.1093/bioinformatics/btab007
