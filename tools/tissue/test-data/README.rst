TISSUE test data
================

The expression fixtures and upstream regression script/tables in this
directory derive from https://github.com/sunericd/TISSUE (``tests/data``
and ``test.py``). See ``LICENSE`` for their license.

The upstream regression test keeps its original settings explicitly:
SpaGE, three folds, ten principal vectors, 15 spatial neighbors, four
gene groups, two cell groups, and alpha 0.23. Its expected tables have
not been changed when updating the wrapper defaults.

Additional fixtures
-------------------

``spatial_custom_key.h5ad`` is ``spatial_sparse.h5ad`` with the coordinates
moved from ``obsm['spatial']`` to ``obsm['vendor_coordinates']``. No
expression values are changed. ``calibration_genes.txt`` lists the first
12 genes of that spatial input, with one repeated name and a blank line
to exercise gene-list parsing.

``calibrated.h5ad`` is generated from the existing ``spatial.h5ad`` and
``reference.h5ad`` with the updated default prediction settings::

    python ../tissue_wrapper.py --spatial spatial.h5ad --reference reference.h5ad \
      --normalize-reference --output calibrated.h5ad --predictions predictions.out.tsv \
      --lower lower.out.tsv --upper upper.out.tsv --uncertainty uncertainty.out.tsv

It then receives an ``obs['condition']`` column with A for the first 64
cells and B for the remaining 64. This annotation only defines the groups
used by downstream tests; it is not a biological condition assignment.

The default bounds and uncertainty tables were updated for the corrected
20-neighbor/automatic-grouping defaults. Prediction values are unchanged.
``knn_predictions.tsv`` uses Harmony-kNN with ten reference neighbors,
three folds, ``reference_normalized.h5ad``, the custom-coordinate input,
and the calibration subset. ``testing.tsv`` uses five imputations, seed
444, and a two-sample t-test of A versus B. ``keep.tsv`` drops the most
uncertain 10 percent of cells independently within A and B (seven cells
per group, using the upstream ceiling rule).

These outputs were generated with TISSUE 1.0.1, harmonypy 0.0.10,
NumPy 1.26.4, SciPy 1.12.0, scikit-learn 1.5.2, and WPCA 0.1.
Numeric regular expressions retain the first three or four decimal
places for numerical-library variation. The unchanged upstream regression
tables retain their original exact comparisons. Every H5AD fixture is
smaller than 1 MB.
