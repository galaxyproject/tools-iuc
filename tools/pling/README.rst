pling
=====

This Galaxy wrapper runs the ``cluster align`` command from `pling`_ to cluster
complete plasmid sequences using sequence-content containment and DCJ-Indel
rearrangement distance.

.. _pling: https://github.com/iqbal-lab-org/pling

User-facing input, parameter, and output instructions are kept in the
``<help>`` section of ``pling_cluster_align.xml``, where they are rendered in
Galaxy.

Implementation notes
^^^^^^^^^^^^^^^^^^^^

The wrapper uses ``write_manifest.py`` to adapt a Galaxy list collection to
Pling's manifest-based command-line interface. The helper:

* stages safe, deterministic FASTA filenames and sorts collection elements by
  identifier;
* detects gzip-compressed FASTA files by content and decompresses them before
  running Pling; and
* validates the optional topology table, requiring one exact ``plasmid`` and
  ``topology`` entry for every collection element and accepting only lowercase
  ``circular`` and ``linear`` values.

The wrapper is limited to complete plasmid sequences. Draft assemblies, short
contigs, and multi-plasmid FASTA files are not supported.

Version note
^^^^^^^^^^^^

The wrapper pins the Conda package to ``pling=3.0.3``. The current package's
``pling --version`` output is ``3.0.2`` because its Python distribution
metadata still reports ``3.0.2``. This is an upstream packaging discrepancy;
the wrapper intentionally follows the Conda package version.

Testing
^^^^^^^

The functional tests cover ordinary and compressed FASTA collection elements,
visualisation outputs, topology handling, sourmash prefiltering, and expected
failures for invalid, duplicate, unknown, and missing topology entries.

License
^^^^^^^

The wrapper is distributed under the MIT license.
