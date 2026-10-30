# Test fixture provenance

Inputs come from the `test_pplace` fixture of nf-core/ampliseq 2.11.0, commit `0473e157ac9a7b1d36845ff9f8fa7ac843d3b234`:
https://github.com/nf-core/ampliseq/blob/0473e157ac9a7b1d36845ff9f8fa7ac843d3b234/conf/test_pplace.config

EPA-ng aligned inputs and placements were retained from the pinned pipeline executed locally in containers. EPA-ng 0.3.8 used the model `GTR+F+I+I+R3` and its default placement-retention settings. GAPPA 0.8.0 produced the expected grafted tree, taxonomy tables, Nexus and PhyloXML outputs. The fully resolved graft fixture was independently generated with the same GAPPA version and `--fully-resolve`.

The SVG regular-expression fixture fixes every line and allows only the generator's wall-clock timestamp to vary. No numerical or graphical difference is ignored. Full execution receipts and sequence-keyed comparisons are retained in the public nekrut/foundry-paper repository under `workflows/ampliseq/maintained-tools/2026-10-29/`.
