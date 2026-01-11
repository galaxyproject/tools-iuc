# SplicingFactory Galaxy tool wrapper

This folder contains a Galaxy wrapper for the R package SplicingFactory providing gene-level splicing diversity calculations.

Files:
- `splicingfactory.xml` - tool definition exposing inputs (transcript matrix and tx2gene mapping), method options, and outputs.
- `splicingfactory_calc.R` - runner script invoked by the tool. Accepts TSV or RDS inputs (tximport or SummarizedExperiment).

Quick test with planemo (requires planemo + conda environment):

```
planemo test tools-iuc/tools/splicingfactory/splicingfactory.xml
```

Notes:
- The tool declares the `splicingfactory` package as a dependency; ensure the package is available to Galaxy via Conda or your environment.
- For TSV inputs the transcript expression matrix must have transcript IDs as row names and a separate `tx2gene` TSV mapping must be provided.
