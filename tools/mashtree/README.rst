Mashtree Galaxy wrapper
=======================

Scope
-----

This wrapper currently targets Mashtree's genome assembly workflow and accepts
FASTA collections. Raw-read, precomputed-sketch, and confidence-estimation
workflows provided by Mashtree are outside the scope of this wrapper.

Implementation notes
--------------------

Input staging
~~~~~~~~~~~~~

Mashtree uses input filenames as sample names. Galaxy datasets therefore need
to be staged with meaningful filenames before Mashtree is executed.

``prepare_inputs.py``:

* treats each collection element as one Mashtree sample;
* derives the Mashtree sample label from the collection element identifier;
* preserves ``.fasta`` and ``.fasta.gz`` suffixes;
* protects identifiers ending in Mashtree-recognized sequence extensions from
  recursive extension stripping;
* detects gzip compression from file content;
* validates basic FASTA structure;
* stages inputs using symbolic links rather than copying sequence data; and
* rejects sample names that cannot be represented uniquely instead of silently
  renaming them.

Mashtree invocation
~~~~~~~~~~~~~~~~~~~

The staged input paths are passed to Mashtree as positional arguments rather
than through ``--file-of-files``.

CPU usage is controlled by Galaxy through ``GALAXY_SLOTS`` and passed to
Mashtree using ``--numcpus``.

Distance matrix ordering
~~~~~~~~~~~~~~~~~~~~~~~~

In testing, Mashtree produced identical distance values for identical
parameters, but the row and column ordering of the output matrix was not
stable between runs. Tests therefore validate matrix contents without relying
on a specific row or column order.

Dependencies
~~~~~~~~~~~~

The wrapper explicitly depends on Python because ``prepare_inputs.py`` is run
as part of the tool command and should not rely on Python being available on
the Galaxy worker host.
