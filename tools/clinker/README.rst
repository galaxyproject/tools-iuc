clinker
-------

This Galaxy wrapper runs `clinker`_ to compare annotated gene clusters or
plasmids and produce a standalone interactive HTML figure.

.. _clinker: https://github.com/gamcil/clinker

Input
^^^^^

Provide an ordered list collection containing at least two annotated GenBank
or GBFF files, one per plasmid. The current Bakta Galaxy wrapper labels its
valid GBFF output as ``tabular``; this wrapper accepts that datatype and parses
each element as GenBank before running clinker. Existing Bakta output can
therefore be connected directly without a second Bakta wrapper or conversion
step.

The identity parameter is a fraction between 0 and 1.
Enable **Keep input collection order** to pass
``--use_file_order``.

Why helper scripts are included
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command is enclosed in XML ``CDATA`` only to keep XML from interpreting
shell characters; CDATA is not an execution language and does not provide a
GenBank parser or structured collection handling. ``prepare_inputs.py`` uses
Biopython to validate every element, require CDS features, preserve collection
order, enforce the two-file minimum, and stage safe filenames for clinker.

Clinker writes input-derived annotation text into a JavaScript JSON object.
``sanitize_html.py`` parses that object as JSON and escapes HTML-sensitive
characters before atomically rewriting the output. Implementing either task as
shell substitutions inside CDATA would be brittle and would not provide the
same validation, fail-closed behavior, or testability.

Output and security
^^^^^^^^^^^^^^^^^^^

The output is a standalone ``figure.html`` file containing the interactive
visualization. Galaxy administrators may need to allow HTML output for this
tool before it can be rendered inline. The wrapper validates the generated
clinker document and escapes ``&``, ``<``, and ``>`` in its embedded JSON, so
GenBank qualifiers cannot close the document's script element.

License and citation
^^^^^^^^^^^^^^^^^^^^

The wrapper is distributed under the MIT license. Please cite:

Gilchrist, C. L. M. and Chooi, Y.-H. (2021). clinker & clustermap.js:
Automatic generation of gene cluster comparison figures.
https://doi.org/10.1093/bioinformatics/btab007
