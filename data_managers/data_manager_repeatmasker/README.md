# RepeatMasker Dfam library data manager

This data manager downloads a [Dfam](https://www.dfam.org/) repeat library in
the partitioned FamDB (v3) format and registers it in the `repeatmasker_famdb`
tool data table for use by the RepeatMasker Galaxy tool (RepeatMasker >= 4.2.4).

## Why this is needed

The RepeatMasker Conda package no longer ships a populated repeat library: the
`Dfam.h5` bundled in the package is a placeholder, and running RepeatMasker
against it fails with `Species "..." is not known to RepeatMasker` because the
library contains no families. The library must be downloaded separately; this
data manager automates that and makes the result available through a data table.

## What it downloads

Dfam distributes its library as a required **root** partition plus optional
**component** partitions, split by curation status (curated/uncurated) and model
type (consensus sequences for the rmblast engine, profile HMMs for nhmmer):

| Component | Engine | Notes |
|-----------|--------|-------|
| Curated consensus | rmblast | Recommended default |
| Uncurated consensus | rmblast | Very large |
| Curated HMM | nhmmer | |
| Uncurated HMM | nhmmer | Very large (100+ partitions) |

The root partition is always downloaded. Each file is verified against its Dfam
md5 checksum, decompressed, and placed in a single directory that RepeatMasker
consumes via `-libdir`.

## Data table

`repeatmasker_famdb` with columns `value, name, version, path`, where `path`
points at the directory holding the FamDB `.h5` partition files.
