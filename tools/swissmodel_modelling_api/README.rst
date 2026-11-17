============================================
SWISS-MODEL Modelling API Galaxy integration
============================================

This is a wrapper for the
`SWISS-MODEL API <https://swissmodel.expasy.org/docs/help#modelling_api>`_.

`SWISS-MODEL <https://swissmodel.expasy.org/>`_ is a web-based integrated
service dedicated to protein structure homology modelling.

With the API wrapper, you are able to build protein models

- fully automated (Automodel mode)
- for a specific alignment provided by the user (Alignment mode)
- based on user provider template coordinates (User template mode)

Models can be created in ModelCIF format and legacy PDB format.

Documentation on how to use the tool inside Galaxy is in the tool XML file.

Installation
------------

Needs Python (>=3.11) and the Python requests package (2.32.5).

Attribution
-----------

When you publish or report results using SWISS-MODEL, please cite the relevant
publications found on https://swissmodel.expasy.org/docs/references (and on the
Galaxy tool page, once installed).

License
-------

License of this Galaxy tool wrapper can be found in `LICENSE <LICENSE>`_.

Terms of use of the SWISS-MODEL server can be found at
https://swissmodel.expasy.org/docs/terms_of_use.

Modelling results of SWISS-MODEL are licensed under the
`CC BY-SA 4.0 Creative Commons Attribution-ShareAlike 4.0 International License <https://creativecommons.org/licenses/by-sa/4.0/legalcode>`_.
