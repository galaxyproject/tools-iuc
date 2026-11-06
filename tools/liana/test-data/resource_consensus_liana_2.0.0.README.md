# LIANA 2.0.0 builtin consensus expectation

resource_consensus_liana_2.0.0.tsv is the complete builtin consensus resource, generated from the released LIANA v2.0.0 data:

- Source: https://github.com/scverse/liana/blob/v2.0.0/src/liana/resource/omni_resource.csv
- Git blob: d004cec380d7a6c971155252c492ad89553e3e0f
- Transformation: the released select_resource("consensus") implementation filters resource == "consensus", selects source_genesymbol/target_genesymbol, renames them ligand/receptor, then writes TSV without an index.
- Implementation: https://github.com/scverse/liana/blob/v2.0.0/src/liana/resource/select_resource.py
- Rows: 4620 (plus header).
- TSV SHA256: 0a089435a1a225b6507679c2079aa30ed234a5c48485278907a3bdcac0ef0d91

The separate resource_consensus.tsv and liana-db/consensus/consensus.tsv fixtures retain the previously cached resource. Cached-resource behavior is tested independently; the builtin expectation deliberately tracks the declared LIANA release.
