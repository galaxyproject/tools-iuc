#!/usr/bin/env python
"""Run mg-toolkit against the EBI Metagenomics v1 API.

mg-toolkit (<=0.10.4) hardcodes its base URL to
``https://www.ebi.ac.uk/metagenomics/api/latest`` in
``mg_toolkit.constants``. That alias now resolves to the v2 API, whose
responses no longer carry the ``meta`` block mg-toolkit parses, so bulk
downloads fail with ``KeyError: 'meta'``.

The v1 API still serves the shape mg-toolkit expects. This wrapper rewrites
the MGnify API URLs to v1 before dispatching to mg-toolkit's normal entry
point, leaving all other behaviour untouched.
"""

import sys

import mg_toolkit.constants as constants
from mg_toolkit import __main__

V1_API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"

# Recompute the MGnify URLs from the v1 base and patch them both in the
# constants module and in the submodules that imported them by value
# (mg_toolkit/__init__ imports the submodules before this runs, so their
# module-level names still point at the old "latest" URLs).
_derived_urls = {
    "API_BASE": V1_API_BASE,
    "MG_SAMPLE_URL": V1_API_BASE + "/samples/{accession}",
    "MG_RUN_URL": V1_API_BASE + "/runs/{accession}",
    "MG_ANALYSES_BASE_URL": V1_API_BASE + "/analyses",
    "MG_ANALYSES_DOWNLOADS_URL": V1_API_BASE + "/analyses/{accession}/downloads",
}
for _name, _value in _derived_urls.items():
    setattr(constants, _name, _value)
for _module_name in ("mg_toolkit.bulk_download", "mg_toolkit.search"):
    _module = sys.modules.get(_module_name)
    if _module is None:
        continue
    for _name, _value in _derived_urls.items():
        if hasattr(_module, _name):
            setattr(_module, _name, _value)

if __name__ == "__main__":
    sys.exit(__main__.main())
