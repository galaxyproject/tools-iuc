Prepare Gemini annotation files and test databases for tool tests
=================================================================

Each version of GEMINI is tied to a particular set of annotation files and
database version.

The ``build-gemini-testdata.sh`` script in this folder should be used to
regenerate the annotation files and the test databases whenever the GEMINI
version required by the tool wrappers gets upgraded.

The script requires a working GEMINI installation at the targeted version and
a folder with GEMINI's original annotation files, and can be executed with::

  sh build-gemini-testdata.sh path/to/gemini/annotation/files
  
It will regenerate the annotation files inside test-data/test-cache/gemini/data
and rebuild the *.db files in test-data.

.. Note::

   If the version of GEMINI that you are upgrading to uses a gemini-config.yaml
   file that is different from the one found in test-data/test-cache you will
   have to upgrade this file manually (make sure you leave the line
   ``annotation_dir: gemini/data`` unchanged in the process).

