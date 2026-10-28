# Synthetic GRN fixtures

`generate_fixtures.py` generates three equivalent 64-cell by 5-gene matrices.
The independent TF_A/TF_B grid drives G1 and G2 respectively, with G3 depending
on both TFs. Values are nonnegative, continuous synthetic expression; these are
software-integration tests, not biological validation or a preprocessing recipe.
The H5AD matrix uses CSR sparse storage. Gene and cell names match the TSV/CSV.

The five Galaxy tests exercise all input readers, known nonempty TF-target
associations, finite positive importance values, absence of self-edges, the
released JSON fitting fields, a nondefault historical stopping configuration
and a no-overlap TF failure with exit code 2. Runtime is excluded from exact
comparison. Tests do not establish biological accuracy or pySCENIC equivalence.

Regenerate with `python generate_fixtures.py` from this directory, using
anndata, numpy, pandas and scipy. HDF5 bytes can vary with their library versions;
the matrix, IDs and text fixtures are deterministic.
