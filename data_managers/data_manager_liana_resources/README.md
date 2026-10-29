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

A single content identity scheme is used for all exported TSVs. Ligand-receptor resources are packaged with LIANA 2.0.0; Metalinks uses an upstream checksum-verified database; HCOP is downloaded from a mutable source without a release identifier in the LIANA API. Hashing the exported TSV distinguishes the actual filtered data even for HCOP and does not use a local download date. The existing four-column table and consumer filename convention are preserved.

Repeated downloads of identical content retain the same identity. Changed content receives a new identity and directory instead of overwriting a previous resource. The display name includes a short digest, while the description records the full source identifier and digest. Existing installed rows, including earlier unversioned cache identifiers, are retained; the migration does not rename or remove them.

The LIANA wrapper uses the registered `path` and appends `<value>.tsv`. Resource-family filtering still uses the stable `value` prefix (`hcop_human_` and `metalinks_`) instead of adding a fifth column.

Development checks
------------------

The standard Galaxy XML tool test checks the consensus resource registration JSON; it does not verify every resource family or deployment into an installed Galaxy data table. The optional unit tests check the downloader's staged files, stable and changed-content identities, and rejection of corrupt staged files. They use the local LIANA stub and pandas, without importing Galaxy internals, and are run separately from Planemo:

```console
python -m unittest discover -s tests -v
```
