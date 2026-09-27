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

Each `path` points to a directory named from the human-readable `name` column and containing `<value>.tsv`. This preserves the layout approved in the original LIANA data-manager review. The LIANA wrapper never reconstructs the directory name: it uses the registered `path` and appends `<value>.tsv`. Resource-family filtering still uses the stable `value` prefix (`hcop_human_` and `metalinks_`) instead of adding a fifth column.
