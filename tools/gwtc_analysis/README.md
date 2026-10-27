# gwtc_analysis admin notes

This file documents deployment requirements for Galaxy administrators.

## Data backends and network requirements

Every analysis that reads catalog products (skymaps, PE files, search-sensitivity injections, LVK data releases)
has a `--data-repo` choice. All backends serve the same files, those of the selected Zenodo release version; a file
a mirror lacks, or a failed transfer, is taken from Zenodo or its public source.

- `zenodo` (default): public Zenodo records and the LIGO DCC; no credentials; outbound HTTPS access required.
- `s3`: the public S3-compatible GWTC bucket used by this tool; no credentials; outbound HTTPS access required.
- `galaxy`: collections staged on the Galaxy instance, else the public usegalaxy.org history of the GWTC products;
  no credentials; outbound HTTPS access required unless all files are staged.

Catalog metadata is looked up through the public GWOSC API (no credentials, outbound HTTPS access required).

## Resources

Most analyses need under 0.5 GB of memory. `search_skymaps` and `parameters_estimation` read full-resolution
sky maps and PE files and peak at about 1.5 GB; `hubble_constant` (bright) at about 1 GB. `rates` and `primordial_black_holes`
read the LVK search-sensitivity injections (~900 MB, downloaded once per job) and peak at about 0.6 GB.

## Not in this wrapper

`galaxy_catalog` and the spectral and dark-siren methods of `hubble_constant` run icarogw, which is not packaged on
conda, and hours-long sampling; they are available from the command line only. `check_catalogs`, a maintenance
check of the package (are new catalogs on GWOSC?), is also left out.

## Dependency

- Conda package: `gwtc_analysis` (version pinned in the wrapper XML).
