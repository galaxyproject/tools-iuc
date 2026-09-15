pling
=====

This Galaxy wrapper runs the ``cluster align`` command from `pling`_ to cluster
complete plasmid sequences using sequence-content containment and DCJ-Indel
rearrangement distance.

.. _pling: https://github.com/iqbal-lab-org/pling

User-facing input, parameter, and output instructions are kept in the
``<help>`` section of ``pling_cluster_align.xml``, where they are rendered in
Galaxy.

Galaxy repository scope
^^^^^^^^^^^^^^^^^^^^^^^

This repository currently provides the ``cluster align`` workflow as a standard
Galaxy tool repository. Additional Pling subcommands can be added later if
there is a need for them.

We considered publishing this as a Galaxy suite, but the current plan is to
provide only ``cluster align``. A suite can be introduced later if additional
Pling subcommands are added.

Implementation notes
^^^^^^^^^^^^^^^^^^^^

The wrapper uses ``write_manifest.py`` to adapt multiple Galaxy FASTA datasets
to Pling's manifest-based command-line interface. The helper:

* stages safe, deterministic FASTA filenames and sorts input datasets by
  identifier;
* detects gzip-compressed FASTA files by content and decompresses them before
  running Pling; and
* validates the optional topology table, requiring one exact ``plasmid`` and
  ``topology`` entry for every selected dataset and accepting only lowercase
  ``circular`` and ``linear`` values. The ``plasmid`` value must include the
  full Galaxy dataset identifier, including a filename extension such as
  ``.fasta`` when it is shown by Galaxy.

The wrapper is limited to complete plasmid sequences. Draft assemblies, short
contigs, and multi-plasmid FASTA files are not supported.

Interactive HTML in Galaxy
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``pling_cluster_align`` tool must be included in Galaxy's HTML sanitisation
allowlist for the interactive plots to work inside Galaxy. The HTML output can
still be downloaded without this, but Galaxy may remove the JavaScript required
for the interactive plots.

Version note
^^^^^^^^^^^^

The wrapper pins the Conda package to ``pling=3.0.3``. The current package's
``pling --version`` output is ``3.0.2`` because its Python distribution
metadata still reports ``3.0.2``. This is an upstream packaging discrepancy;
the wrapper intentionally follows the Conda package version.

Testing
^^^^^^^

The functional tests cover ordinary and compressed FASTA datasets,
visualisation outputs, topology handling, sourmash prefiltering, and expected
failures for invalid, duplicate, unknown, and missing topology entries.

License
^^^^^^^

The wrapper is distributed under the MIT license.
