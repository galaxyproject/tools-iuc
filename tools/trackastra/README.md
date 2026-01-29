# Trackastra Galaxy Tool

## Overview

**Trackastra** is a deep learning-based tool for tracking cell instances in time-lapse microscopy images. It combines:
- **Cellpose**: For automatic cell segmentation
- **Trackastra**: For transformer-based cell tracking across time

This Galaxy wrapper enables easy access to cell tracking workflows without requiring command-line expertise.

## Requirements

### Input Data

- **Zarr dataset**: Time-series images in zarr format (NGFF-compliant)
  - Must have dimensions: `t` (time), `z` (depth), `y` (rows), `x` (columns)
  - Optional: extra dimensions (e.g., channels, which can be selected via coordinates)
  - Data type: typically uint16 (16-bit unsigned integer)

### Software Dependencies

Managed via `pixi.toml`:
- `trackastra >= 0.3.2`
- `cellpose < 4.0`
- `ngff_zarr >= 0.20.0`
- `s3fs >= 2026.1.0` (for remote data access)
- `napari` (for visualization)

## Usage Modes

### 1. Segment and Track (Recommended)

Automatically segments cells using Cellpose v3, then performs tracking.

```bash
python trackastra_wrapper.py segment_and_track \
    --zarr_path https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0101A/13457537.zarr/0 \
    --scale_level 0 \
    --channel_coords 0 \
    --downscale_x 1.0 \
    --downscale_y 1.0 \
    --downscale_z 1.0 \
    --start_tp 0 \
    --end_tp -1 \
    --segmentation_model cyto3 \
    --tracking_model ctc
```

The tool creates two output pickle files with hardcoded names:
- `trackastra_trackgraph.pickle.dat` - Native Trackastra tracking object
- `napari_tracks.pickle.dat` - Napari-compatible tracking data

### 2. Track Only

Assumes pre-existing segmentation in a separate channel/dimension.

```bash
python trackastra_wrapper.py track \
    --zarr_path https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0101A/13457537.zarr/0 \
    --raw_channel_coords 0 \
    --seg_channel_coords 1 \
    --downscale_x 1.0 \
    --downscale_y 1.0 \
    --downscale_z 1.0 \
    --start_tp 0 \
    --end_tp -1 \
    --tracking_model ctc
```

## Parameters

### Common Parameters (Both Modes)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `zarr_path` | Required | URL (https://..., s3://...) or local path to OME-Zarr dataset |
| `scale_level` | 0 | Pyramid level in zarr (0 = finest/best resolution) |
| `downscale_x`, `downscale_y`, `downscale_z` | 1.0 | Spatial downscaling factors (>1 reduces resolution for speed) |
| `start_tp` | 0 | First time frame to process (0-indexed) |
| `end_tp` | -1 | Last time frame to process (-1 = all frames) |
| `tracking_model` | ctc | Trackastra tracking model identifier |
| `channel_coords` / `raw_channel_coords` | 0 | Non-tzyx dimension coordinates (space or comma-separated) |

### Segment and Track Mode Only

| Parameter | Default | Description |
|-----------|---------|-------------|
| `segmentation_model` | cyto3 | Cellpose v3 model: `cyto3`, `cyto2`, or `nuclei` |
| `channel_coords` | 0 | Coordinate to select raw image channel in multi-dimensional zarr |

### Track Only Mode Only

| Parameter | Default | Description |
|-----------|---------|-------------|
| `raw_channel_coords` | 0 | Coordinate(s) to select raw image channel |
| `seg_channel_coords` | 0 | Coordinate(s) to select segmentation channel |

### Multi-Dimensional Zarr Navigation

For zarr datasets with extra dimensions beyond `tzyx`:
- Extra dimensions are automatically moved to the front
- Provide integer coordinates for each extra dimension
- Example: 6D data `(t, view, domain, z, y, x)` needs coordinates like `"0 1"` to select view=0, domain=1
- Use space or comma separation: `"0 1"` or `"0,1"`

## Output Files

The tool generates **two pickle files** with fixed names:

1. **trackastra_trackgraph.pickle.dat**: Native Trackastra tracking object
   - Contains the complete tracking graph with all structural information
   - Can be loaded and analyzed programmatically with pickle

2. **napari_tracks.pickle.dat**: Napari-compatible tracking data
   - Pre-formatted for visualization in Napari
   - Load in Napari with: `viewer.add_tracks(pickle.load(f))`

Both files are created in the working directory and automatically copied to the Galaxy outputs.

## Remote Data Access

The tool supports various zarr data sources:

```bash
# Activate environment
pixi shell

# Run segment_and_track
python trackastra_wrapper.py segment_and_track \
    --zarr_path test-data/sample_timelapse.zarr \
    --end_tp 3 \
    --output_tracks test_output.csv

# View results
cat test_output.csv
```

## Implementation Notes

### Tool Wrapper Structure

```
tools/trackastra/
├── trackastra.xml           # Galaxy tool definition (XML)
├── trackastra_wrapper.py    # Main Python CLI wrapper
├── pixi.toml                # Pixi environment dependencies
├── .shed.yml                # Tool Shed metadata
├── test-data/               # Test inputs directory (initially empty)
└── README.md                # This file
```

## Galaxy Integration

### Command Line Interface

The tool provides two subcommands for Galaxy:

**Segment and Track:**
```
python trackastra_wrapper.py segment_and_track --zarr_path ... --channel_coords ... --segmentation_model cyto3 ...
```

**Track Only:**
```
python trackastra_wrapper.py track --zarr_path ... --raw_channel_coords ... --seg_channel_coords ... ...
```

All coordinates and numerical parameters are validated before execution. Errors print to stderr with appropriate exit codes for Galaxy error detection.

## Online Zarr Datasets

Ready-to-use test datasets from IDR (Image Data Resource):

```
https://uk1s3.embassy.ebi.ac.uk/idr/zarr/v0.4/idr0101A/13457537.zarr/0
```

This dataset can be used directly as zarr_path for testing without downloading.

## Advanced: Key Parameters Explained

**Downscaling**: 
- Useful for large images (500+ pixels per dimension)
- Example: `--downscale_x 2.0` processes at half resolution, faster but may reduce tracking precision
- Set to 1.0 for full resolution (default)
- Consider downscaling when processing large 3D datasets or memory-limited systems

**Time Point Range**:
- Helpful for testing parameters on a subset
- `--start_tp 0 --end_tp 5` processes frames 0-5 only
- Use `-1` for end_tp to include all frames
- Useful for validating parameters before processing entire timeseries

**Pyramid Levels**:
- Zarr datasets often have multi-resolution pyramids
- Level 0 = finest resolution (slowest, most accurate)
- Higher levels = downsampled versions (faster)
- Default level 0 is recommended

**Channel/Dimension Coordinates**:
- If zarr has extra dimensions beyond tzyx (e.g., tchannel,z,y,x), use coordinate indices
- Example: for tchannel,z,y,x format use `--channel_coords 0` to select first channel
- Multiple coordinates use space or comma separation: `0 1` or `0,1`

## Troubleshooting

### Common Issues and Solutions

**"Scale index negative or larger than available resolutions"**
- The requested pyramid level doesn't exist in the zarr
- Solution: Use `--scale_level 0` (default, safest option)

**Memory errors on large datasets**
- The downscaled data is still too large for available memory
- Solutions:
  - Increase `--downscale_x/y/z` factors further
  - Process fewer time points with `--start_tp` and `--end_tp`
  - Use a higher `--scale_level` (lower resolution)

**Segmentation quality issues**
- Segmentation is too aggressive or too lenient
- Try different `--segmentation_model` options:
  - `cyto3`: General cells (recommended, default)
  - `cyto2`: Alternative model for comparison
  - `nuclei`: Optimized for nuclear/DAPI staining
- Check image contrast and brightness

**Tracking breaks or incomplete tracks**
- Cells moving too fast between frames may be missed
- Solutions:
  - Reduce `--downscale_*` factors for finer tracking
  - Check image quality and contrast
  - Reduce frame intervals if possible (i.e., use more timepoints)

### Debugging

Enable verbose output:
```bash
python -u trackastra_wrapper.py segment_and_track \
    --zarr_path data.zarr \
    --output_tracks tracks.csv 2>&1 | tee debug.log
```

## Citation

If you use Trackastra in published research, please cite:
- Trackastra: [GitHub](https://github.com/QuantumAstronomy/trackastra)
- Cellpose: Stringer, C., Wang, T., Michaelos, M., & Pachitariu, M. (2021). Cellpose: a generalist algorithm for cellular segmentation. *Nature Methods*, 18(1), 100-106.

## References

- [NGFF Zarr Spec](https://ngff.openmicroscopy.org/)
- [Trackastra Documentation](https://github.com/QuantumAstronomy/trackastra)
- [Cellpose Documentation](https://cellpose.readthedocs.io/)
- [Galaxy Tool Development](https://docs.galaxyproject.org/en/master/dev/schema.html)
