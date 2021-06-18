# Data preparation

### Source

All GISAID data is downloaded and run through [`grapevine`](https://github.com/cov-ert/grapevine) which excludes records without proper dates, removes duplicate sequences (taking the earliest sample of the duplicates), omits some sequences with known issues, filters by length and coverage, and trims the sequences to CDS.

It also aligns the sequences using `mafft` and builds an ML tree using `iqtree`. A lineages is assigned to each sequence using `pangolin` with the previous data release.

### Lineage Curation

The phylogeny is annotated with lineage and then in `FigTree` the lineages are manually curated, drawing together a number of pieces of information including monophyly in the ML phylogeny (generally a bootstrap > 70 is required) and epidemiological data such as country and travel history. Any changes to lineage definitions and new lineages are documented during this process.

- The lineage may have been defined earlier in the outbreak and with added sequence data, there is less support for that lineage. In these cases the associated epidemiological metadata is examined and the lineage may be refined or even dropped entirely. The lineage number will not be 'recycled', but the members will get reassigned the parent lineage designation.
- The lineage may have very clear epidemiological support and ambiguities or homoplasies in the sequences/ tree could contribute to low bootstrap values. In these cases, if the support is strong, the lineages are called. Recall rates for these lingeages within `pangolin` may be lower however.


