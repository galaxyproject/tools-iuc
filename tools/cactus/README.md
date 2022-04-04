# Cactus

This suite contains a wrapper for the [Cactus](github.com/comparativeGenomicsToolkit/cactus) whole-genome aligner and another wrapper for exporting the data it generates.

## Requirements

The developers provide an official Docker container.
However, the Docker container won't run unless `/etc/passwd` is mounted ([see this issue](https://github.com/ComparativeGenomicsToolkit/cactus/issues/677)).
Because of this, it is recommended to run Cactus using Singularity.

Cactus uses a lot of RAM. We have tested it on Galaxy using 24 GB of RAM for a
progressive alignment of three chromosome-level, 150 megabase genomes.