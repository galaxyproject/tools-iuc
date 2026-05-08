# ProSeg Galaxy Tool Integration

This repository contains the Galaxy tool wrapper for **ProSeg** (Probabilistic Cell Segmentation), developed as part of the **Spatial2Galaxy / EISTA** workflow integration.

## Current Progress
- [x] **XML Wrapper:** Fully developed with appropriate input/output parameters.
- [x] **Linting:** Passed `planemo lint` (all requirements met).
- [x] **Functional Testing:** Passed `planemo test` using a custom 3D transcript dataset.
- [x] **Version Control:** Verified and pushed to GitHub.

## Technical Details
- **Tool ID:** `proseg`
- **Primary Outputs:** Zarr-based segmentation results.
- **Dependencies:** ProSeg 3.1.1 (via Conda).

## Next Steps
- Validate integration within the Freiburg EISTA environment.
- Test with real-world spatial transcriptomics datasets (Xenium/Molecular Cartography).

**Contributor:** Firas Manasrah (Charité BIH)
