# gwtc_analysis admin notes

This file documents deployment requirements for Galaxy administrators.

## Data backends and network requirements

- `galaxy`: reads data available on the Galaxy instance; no outbound network access required.
- `s3`: reads from the public S3-compatible GWTC bucket used by this tool; no credentials required; outbound HTTPS access required.
- `zenodo`: downloads from public Zenodo records; no credentials required; outbound HTTPS access required.

For `catalog_statistics` and `event_selection`, metadata lookup via the public GWOSC API may be performed; no credentials are required, but outbound HTTPS access is still required.

## Dependency

- Conda package: `gwtc_analysis` (version pinned in the wrapper XML).
