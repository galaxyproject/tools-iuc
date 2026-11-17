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
  order, and stages safe filenames for clinker.
* ``sanitize_html.py`` parses clinker's embedded JSON and escapes
  HTML-sensitive characters before rewriting the output, preventing
  input-derived annotation text from closing the document's script element.

HTML output
^^^^^^^^^^^

Galaxy administrators may need to allow HTML output for this tool before the
generated visualization can be rendered inline.

Testing
^^^^^^^

The functional tests cover ordered collections, generic annotated records,
invalid input, the one-file failure case, and a regression test for
HTML/JavaScript injection in annotation text.

License
^^^^^^^

The wrapper is distributed under the MIT license.
