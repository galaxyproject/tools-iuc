JBrowse in Galaxy
=================

    JBrowse is a fast, embeddable genome browser built completely with
    JavaScript and HTML5

Thus, it makes an ideal fit with Galaxy, especially for use as a
workflow summary. E.g. annotate a genome, then visualise all of the
associated datasets as an interactive HTML page. This tool MUST be whitelisted
(or ``sanitize_all_html=False`` in galaxy.ini) to function correctly.

Installation
============

It is recommended to install this wrapper via the Galaxy Tool Shed.

Running Locally
===============

The Galaxy tool interface writes out a xml file which is then used to generate
the visualizations. An example used during development/testing can be seen in
`test-data/*/test.xml`. The format is in no way rigorously defined and is
likely to change at any time! Beware. ;)

History
=======

- 1.16.5+galaxy0

    - UPDATED to JBrowse 1.16.5

- 1.16.4+galaxy0

    - UPDATED to JBrowse 1.16.4
    - ADDED filter too big metadata
    - CHANGED default value for topLevelFeatures (gene subfeatures are now inferred) and style.className (feature style was fixed)

- 1.16.2+galaxy0

    - UPDATED to JBrowse 1.16.2
    - ADDED support for NeatHTMLFeatures and NeatCanvasFeatures track types

- 1.16.1+galaxy0

    - UPDATED to JBrowse 1.16.1
    - ADDED support for MultiBigWig plugin
    - ADDED support for tabix indexing of fasta and gff
    - ADDED support for REST and SPARQL endpoints
    - ADDED option to change chunk size for BAM tracks
    - FIXED loading of VCF files. They were gzipped and the URLs were incorrect
    - FIXED metadata on tracks types other than GFF+HTML
    - FIXED infrastructure URL parsing (and embedding in links) for some tracks
    - REMOVED support for selecting multiple genomes as input due to tracking of track metadata
    - REMOVED support for themes as JBrowse no longer allow runtime loading of plugins

- 0.7 Support for plugins (currently GC Content, Bookmarks, ComboTrackSelector),
  track metadata
- 0.5.2 Support for CanvasFeatures options.
- 0.5.1 Support for contextual menus. Conda tests.
- 0.5 Update existing instances on disk. Index names. Support HTML tracks
  instead of Canvas. Support default tracks. General JBrowse optinos
- 0.4 Support for dataset collections and customisation of tracks including
  labelling, colours, styling. Added support for genetic code selection.
  Fixed package installation recipe issues.
- 0.3 Added support for BigWig, etc.
- 0.2 Added support for BAM, Blast, VCF.
- 0.1 Initial public release.

Wrapper License (MIT/BSD Style)
===============================

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the names
of the contributors or copyright holders not be used in advertising or
publicity pertaining to distribution of the software without specific
prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
