#!/usr/bin/env python3
"""Convert Cell Tracking Challenge results to Napari-friendly outputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

import tifffile
from skimage.measure import regionprops

MASK_PATTERN = re.compile(r"^man_track(?P<frame>[0-9]+)\.tif$")


@dataclass(frozen=True)
class TrackRecord:
    """One row from a CTC man_track.txt table."""

    track_id: int
    start_frame: int
    end_frame: int
    parent_track_id: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert CTC masks and man_track.txt to a Napari tracks table, "
            "a lineage graph, and a portable CTC archive."
        )
    )
    parser.add_argument("--track-table", required=True, type=Path)
    parser.add_argument("--masks-dir", required=True, type=Path)
    parser.add_argument("--tracks-out", required=True, type=Path)
    parser.add_argument("--graph-out", required=True, type=Path)
    parser.add_argument("--archive-out", required=True, type=Path)
    return parser.parse_args()


def read_track_table(path: Path) -> list[TrackRecord]:
    if not path.is_file():
        raise ValueError(f"CTC track table does not exist: {path}")

    records: list[TrackRecord] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(
                f"CTC track-table line {line_number} must contain exactly "
                f"four columns; found {len(fields)}."
            )
        try:
            values = [int(field) for field in fields]
        except ValueError as exc:
            raise ValueError(
                f"CTC track-table line {line_number} contains a non-integer value."
            ) from exc
        records.append(TrackRecord(*values))

    track_ids = [record.track_id for record in records]
    if len(track_ids) != len(set(track_ids)):
        raise ValueError("CTC track IDs must be unique in man_track.txt.")
    if any(record.track_id <= 0 for record in records):
        raise ValueError("CTC track IDs must be positive integers.")
    if any(
        record.start_frame < 0 or record.end_frame < 0
        for record in records
    ):
        raise ValueError("CTC frame numbers must be non-negative integers.")
    if any(record.start_frame > record.end_frame for record in records):
        raise ValueError("A CTC track cannot end before it starts.")

    track_id_set = set(track_ids)
    parent_ids = {
        record.parent_track_id
        for record in records
        if record.parent_track_id != 0
    }
    missing_parents = parent_ids - track_id_set
    if missing_parents:
        raise ValueError(
            "CTC parent track IDs are missing from the track table: "
            + ", ".join(str(value) for value in sorted(missing_parents))
        )
    if any(
        record.track_id == record.parent_track_id
        for record in records
    ):
        raise ValueError("A CTC track cannot be its own parent.")

    records_by_id = {record.track_id: record for record in records}
    for record in records:
        if record.parent_track_id == 0:
            continue
        parent = records_by_id[record.parent_track_id]
        if parent.end_frame >= record.start_frame:
            raise ValueError(
                f"Parent track {parent.track_id} must end before child track "
                f"{record.track_id} starts."
            )

    return sorted(records, key=lambda record: record.track_id)


def find_masks(directory: Path) -> list[tuple[int, Path]]:
    if not directory.is_dir():
        raise ValueError(f"CTC masks directory does not exist: {directory}")

    masks: list[tuple[int, Path]] = []
    for path in directory.iterdir():
        if not path.is_file():
            continue
        match = MASK_PATTERN.fullmatch(path.name)
        if match:
            masks.append((int(match.group("frame")), path))

    if not masks:
        raise ValueError(f"No CTC mask TIFF files were found in {directory}.")

    masks.sort(key=lambda item: item[0])
    expected_frames = list(range(len(masks)))
    actual_frames = [frame for frame, _ in masks]
    if actual_frames != expected_frames:
        raise ValueError(
            "CTC mask frames must be consecutive and start at 0; found "
            + ", ".join(str(frame) for frame in actual_frames)
        )
    return masks


def build_napari_outputs(
    table: list[TrackRecord],
    masks: list[tuple[int, Path]],
) -> tuple[list[str], list[list[float | int]], dict[int, list[int]]]:
    rows: list[list[float | int]] = []
    spatial_ndim: int | None = None
    expected_shape: tuple[int, ...] | None = None
    expected_dtype = None

    for frame, path in masks:
        mask = tifffile.imread(path)
        if mask.ndim not in (2, 3):
            raise ValueError(
                f"Expected a 2D or 3D CTC mask in {path.name}; found {mask.shape}."
            )
        if spatial_ndim is None:
            spatial_ndim = mask.ndim
            expected_shape = mask.shape
            expected_dtype = mask.dtype
        elif mask.ndim != spatial_ndim:
            raise ValueError("All CTC masks must have the same dimensionality.")
        elif mask.shape != expected_shape:
            raise ValueError("All CTC masks must have the same spatial shape.")
        elif mask.dtype != expected_dtype:
            raise ValueError("All CTC masks must have the same data type.")

        active_ids = {
            record.track_id
            for record in table
            if record.start_frame <= frame <= record.end_frame
        }
        regions = list(regionprops(mask))
        observed_ids = {int(region.label) for region in regions}
        if observed_ids != active_ids:
            missing = active_ids - observed_ids
            unexpected = observed_ids - active_ids
            details = []
            if missing:
                details.append(
                    "missing labels " + ", ".join(str(value) for value in sorted(missing))
                )
            if unexpected:
                details.append(
                    "unexpected labels "
                    + ", ".join(str(value) for value in sorted(unexpected))
                )
            raise ValueError(
                f"CTC table and mask {path.name} disagree: " + "; ".join(details)
            )

        for region in regions:
            rows.append(
                [int(region.label), frame, *[float(value) for value in region.centroid]]
            )

    if spatial_ndim == 2:
        header = ["track_id", "t", "y", "x"]
    elif spatial_ndim == 3:
        header = ["track_id", "t", "z", "y", "x"]
    else:
        raise RuntimeError("Could not determine CTC mask dimensionality.")

    rows.sort(key=lambda row: (int(row[0]), int(row[1])))
    graph = {
        record.track_id: [record.parent_track_id]
        for record in table
        if record.parent_track_id != 0
    }
    return header, rows, graph


def write_tracks(path: Path, header: list[str], rows: list[list[float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_graph(path: Path, graph: dict[int, list[int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(graph, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_archive(
    path: Path,
    track_table: Path,
    masks: list[tuple[int, Path]],
    tracks_out: Path,
    graph_out: Path,
) -> None:
    readme = """Trackastra visualization bundle

tracked_ctc/
  man_track.txt          Cell Tracking Challenge lineage table
  man_trackNNNN.tif      Relabelled CTC masks
napari_tracks.csv        Napari Tracks vertices: track_id, t, (z), y, x
napari_graph.json        Child track ID to parent track ID list

The tracked_ctc directory can be opened with napari-trackastra. The CSV and
JSON files can also be loaded programmatically::

    import json
    import pandas as pd
    import napari

    tracks = pd.read_csv("napari_tracks.csv").to_numpy()
    with open("napari_graph.json") as handle:
        raw_graph = json.load(handle)
    graph = {int(child): [int(parent) for parent in parents]
             for child, parents in raw_graph.items()}

    viewer = napari.Viewer()
    viewer.add_tracks(tracks, graph=graph)
    napari.run()
"""

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="trackastra_export_") as temporary:
        root = Path(temporary)
        ctc_dir = root / "tracked_ctc"
        ctc_dir.mkdir()
        shutil.copyfile(track_table, ctc_dir / "man_track.txt")
        for _, mask in masks:
            shutil.copyfile(mask, ctc_dir / mask.name)
        shutil.copyfile(tracks_out, root / "napari_tracks.csv")
        shutil.copyfile(graph_out, root / "napari_graph.json")
        (root / "README.txt").write_text(readme, encoding="utf-8")

        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for file_path in sorted(root.rglob("*")):
                if file_path.is_file():
                    archive.write(file_path, file_path.relative_to(root))


def main() -> None:
    args = parse_args()
    table = read_track_table(args.track_table)
    masks = find_masks(args.masks_dir)
    header, rows, graph = build_napari_outputs(table, masks)
    write_tracks(args.tracks_out, header, rows)
    write_graph(args.graph_out, graph)
    write_archive(
        args.archive_out,
        args.track_table,
        masks,
        args.tracks_out,
        args.graph_out,
    )


if __name__ == "__main__":
    main()
