This data manager downloads and installs Savont reference databases so they
can be used by the **Savont classify** and **Savont SINTAX** tools
(`tools/savont <https://github.com/galaxyproject/tools-iuc/tree/main/tools/savont>`_).

The available databases (see the
`Savont documentation <https://github.com/bluenote-1577/savont>`_) are:

- **emu-1**: EMU 16S rRNA reference database.
- **silva-138.2**: SILVA 16S rRNA reference database (v138.2).
- **greengenes2-2024.09**: GreenGenes2 16S rRNA reference database (2024.09).
- **unite-10.0**: UNITE ITS reference database (v10.0).

Databases can be several gigabytes in size; make sure the Galaxy installation
has enough free disk space under ``data/data_managers`` before running this
data manager.

Tests install a small mock database (a re-labelled subset of the Savont
reference ASVs) instead of downloading a full database.
