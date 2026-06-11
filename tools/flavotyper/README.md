# FlavoTyper Galaxy tool

Galaxy wrapper for [FlavoTyper](https://forge.inrae.fr/eric.duchaud/flavotyper) —
in silico serotyping of *Flavobacterium psychrophilum* from genome assemblies.

FlavoTyper detects specific biomarkers in a genome assembly with BLASTN, verifies the
species with fastANI, and assigns the serotype. Tool dependencies are resolved through
the [`flavotyper` Bioconda package](https://anaconda.org/bioconda/flavotyper).

## Contents

| File | Purpose |
|------|---------|
| `flavotyper.xml` | Tool wrapper |
| `macros.xml` | Version token, Bioconda requirement, citation |
| `test-data/` | Example genome used by the functional tests |

## Development

```bash
planemo lint flavotyper.xml
planemo test flavotyper.xml
```
