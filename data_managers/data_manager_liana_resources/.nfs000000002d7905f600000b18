LIANA Resources Data Manager
============================

This Galaxy data manager downloads and caches resources used by the LIANA+ 1.8.1 wrappers.

Managed resource families
-------------------------

- LIANA ligand-receptor resources returned by `li.resource.show_resources()`
- HCOP human-to-target-organism ortholog maps
- Metalinks metabolite-protein interaction presets

The `liana_resources` data table intentionally retains the existing four-column schema:

1. `value`
2. `name`
3. `description`
4. `path`

Each `path` points to a directory containing `<value>.tsv`. The LIANA wrapper filters entries by the `value` prefix (`hcop_human_` and `metalinks_`) instead of adding a fifth resource-type column.
