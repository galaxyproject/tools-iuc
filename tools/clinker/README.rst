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

Producing input with Bakta
~~~~~~~~~~~~~~~~~~~~~~~~~~

When using the current Bakta Galaxy wrapper, keep CDS annotation enabled (do
not select ``--skip-cds``) and select **Annotations and sequences in GenBank
format** (``file_gbff``) in **Selection of the output files**. Bakta exposes
this output as ``annotation_gbff`` and Galaxy may label it ``tabular``; that is
expected and is accepted by this wrapper. When Bakta is run on a FASTA
collection, use the resulting ``annotation_gbff`` list collection as the
clinker input. When Bakta is run one dataset at a time, collect at least two
GBFF outputs into an ordered list before running clinker. The other Bakta
settings (database, organism metadata, complete-replicon mode, and so on) can
be chosen according to the annotation task; they are not clinker-specific.

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

Saving and restoring interactive edits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The rendered figure provides **Save Data** and **Save SVG** controls. **Save
Data** downloads ``clinker.json``, containing the current plot data and
configuration (for example renamed labels, colours, hidden links, and layout
options). **Save SVG** downloads a static vector image of the current figure.
Galaxy's dataset download remains the original generated HTML; browser edits
are not written back into that dataset.

To preserve an editable result, keep the original HTML together with the
saved ``clinker.json``. Open the HTML locally or in Galaxy, use **Load data**,
and select the JSON file. Uploading the JSON to the Galaxy history is useful
for provenance and sharing, but the browser file chooser still needs a local
copy when loading it into the visualization.

License and citation
^^^^^^^^^^^^^^^^^^^^

The wrapper is distributed under the MIT license. Please cite:

Gilchrist, C. L. M. and Chooi, Y.-H. (2021). clinker & clustermap.js:
Automatic generation of gene cluster comparison figures.
https://doi.org/10.1093/bioinformatics/btab007
