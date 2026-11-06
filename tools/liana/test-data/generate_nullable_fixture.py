#!/usr/bin/env python
"""Generate the nullable-string regression input without changing expression.

Source: https://zenodo.org/records/18388645/files/pbmc_small.h5ad
Input SHA256:
f2b88dfbb0fcb3674e36c1dcf3573e8f8d1a8f26c26da66684e33a1a134b1a6b
Generated with AnnData 0.13.0 and pandas 3.1.0, with gzip compression.
Run with explicit --input and --output paths in an AnnData >=0.11 environment.
The original matrices, raw expression and annotations are retained. Only the
observation index dtype/name and a new nullable_note column are changed.
"""

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def generate_fixture(input_path, output_path):
    adata = ad.read_h5ad(input_path)
    string_dtype = pd.StringDtype(storage="python")
    adata.obs.index = pd.Index(
        adata.obs.index.to_numpy(dtype=object), dtype=string_dtype, name=None
    )
    adata.obs["nullable_note"] = pd.array(
        [pd.NA if i == 0 else f"fixture_note_{i}" for i in range(adata.n_obs)],
        dtype=string_dtype,
    )
    with ad.settings.override(allow_write_nullable_strings=True):
        adata.write_h5ad(
            output_path,
            compression="gzip",
            convert_strings_to_categoricals=False,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate_fixture(args.input, args.output)


if __name__ == "__main__":
    main()
