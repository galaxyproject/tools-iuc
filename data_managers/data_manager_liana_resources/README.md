LIANA Resources Data Manager
============================

This Galaxy data manager downloads and caches resources used by the LIANA+ 2.0.0 wrappers.

Managed resource families
-------------------------

- LIANA ligand-receptor resources returned by `li.rs.show_resources()`
- HCOP human-to-target-organism ortholog maps
- Metalinks metabolite-protein interaction presets

The `liana_resources` data table intentionally retains the existing four-column schema:

1. `value`
2. `name`
3. `description`
4. `path`

Downloaded files are staged under the output dataset's `extra_files_path`, so Galaxy can move them into its managed data directory. Each new `value` combines the source resource identifier with the full SHA256 digest of its TSV content, for example `consensus_sha256_<digest>`. Its installed `path` is `liana-db/<value>` and contains `<value>.tsv`.

Repeated downloads of identical content retain the same identity. Changed content receives a new identity and directory instead of overwriting a previous resource. The display name includes a short digest, while the description records the full source identifier and digest. Existing installed rows, including earlier unversioned cache identifiers, are retained; the migration does not rename or remove them.

The LIANA wrapper uses the registered `path` and appends `<value>.tsv`. Resource-family filtering still uses the stable `value` prefix (`hcop_human_` and `metalinks_`) instead of adding a fifth column.

Development checks
------------------

The Galaxy tool test checks the resource registration JSON. The additional move-aware tests exercise Galaxy's actual data-manager file processor, verify the registered file exists after moving, preserve an earlier installation, and check stable and changed-content identities. Run them in a development environment containing LIANA's dependencies and `galaxy-tool-util` (as provided by Planemo):

```console
python -m unittest discover -s tests -v
```
