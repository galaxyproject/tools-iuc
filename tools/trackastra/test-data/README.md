# Trackastra test data

- `image_000*.tif` and `mask_000*.tif`: three-frame linear tracking fixture.
- `division_images_stack.tif` and `division_masks_stack.tif`: multi-page TIFF
  fixture containing one parent and two daughter objects.
- `expected_linear_000*.tif`: exact tracked labels for the linear fixture.
- `expected_division_000*.tif`: exact CTC track IDs for the division fixture.
- `tiny_model/`: deterministic Trackastra model used only by tests.
- `trackastra_models.loc`: test Tool Data Table entry for `tiny_model/`.

The edge and lineage tables are asserted inline in `trackastra.xml`; no
separate expected CSV/TXT fixtures are required.
