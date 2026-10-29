#!/usr/bin/env python
"""Generate constant spatial curves with analytically known AUC scores.

For radii 0, 1, 2, g=2 has mean log2 score 1 and g=0.5 has score -1.
Both curves peak at the first radius because their absolute values tie.
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def generate_fixture(output_path):
    rows = []
    for source, target, ligand, receptor, value in [
        ("A", "B", "L1", "R1", 2.0),
        ("B", "A", "L2", "R2", 0.5),
    ]:
        for radius in [0.0, 1.0, 2.0]:
            rows.append({
                "source": source,
                "target": target,
                "ligand_complex": ligand,
                "receptor_complex": receptor,
                "interaction": f"{ligand}^{receptor}",
                "radius": radius,
                "g": value,
                "g_expr": value,
                "g_pcf": 1.0,
            })
    adata = ad.AnnData(np.ones((1, 1)))
    adata.uns["lric"] = pd.DataFrame(rows)
    with ad.settings.override(allow_write_nullable_strings=True):
        adata.write_h5ad(output_path, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    generate_fixture(parser.parse_args().output)
