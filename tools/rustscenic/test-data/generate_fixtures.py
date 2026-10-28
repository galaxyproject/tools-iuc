"""Create the small synthetic, format-equivalent GRN integration fixtures.

These matrices test tool plumbing, not biological accuracy. No random numbers
or external data are used. Run with anndata, numpy, pandas and scipy installed.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).parent


def main():
    """Write equivalent cells-by-genes TSV, CSV and sparse AnnData inputs."""
    cells = np.arange(64)
    tf_a = (cells % 8).astype(np.float32) / 4
    tf_b = (cells // 8).astype(np.float32) / 4
    values = np.column_stack(
        [tf_a, tf_b, 1 + 2 * tf_a, 1 + 3 * tf_b, (tf_a + tf_b) / 2]
    ).astype(np.float32)
    frame = pd.DataFrame(
        values,
        index=pd.Index([f"cell_{i:02d}" for i in cells], name="cell"),
        columns=["TF_A", "TF_B", "G1", "G2", "G3"],
    )
    frame.to_csv(ROOT / "expression.tsv", sep="\t")
    frame.to_csv(ROOT / "expression.csv")
    adata = ad.AnnData(
        X=sparse.csr_matrix(values),
        obs=pd.DataFrame(index=frame.index),
        var=pd.DataFrame(index=frame.columns),
    )
    adata.write_h5ad(ROOT / "expression.h5ad")
    (ROOT / "tfs.txt").write_text("# Candidate transcription factors\nTF_A\n\nTF_B\n")
    (ROOT / "missing_tfs.txt").write_text("ABSENT_TF\n")


if __name__ == "__main__":
    main()
