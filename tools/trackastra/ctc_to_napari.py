#!/usr/bin/env python3
"""Convert Cell Tracking Challenge results to Napari-friendly outputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from skimage.measure import regionprops
from trackastra.tracking import ctc_to_napari_tracks

MASK_PATTERN = re.compile(r"^man_track(?P<frame>[0-9]+)\.tiff?$")
GALAXY_TIFF_SUFFIX = ".tiff"


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
    parser.add_argument("--stack-out", required=True, type=Path)
    parser.add_argument("--tracks-out", required=True, type=Path)
    parser.add_argument("--graph-out", required=True, type=Path)
    parser.add_argument("--ctc-dir", required=True, type=Path)
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
) -> tuple[
    list[str],
    list[list[float | int]],
    dict[int, list[int]],
    "np.ndarray",
]:
    frames: list["np.ndarray"] = []
    reference: "np.ndarray" | None = None

    for frame, path in masks:
        mask = tifffile.imread(path)
        if mask.ndim not in (2, 3):
            raise ValueError(
                f"Expected a 2D or 3D CTC mask in {path.name}; found {mask.shape}."
            )
        if reference is None:
            reference = mask
        elif mask.shape != reference.shape or mask.dtype != reference.dtype:
            raise ValueError(
                "All CTC masks must have the same shape and data type."
            )

        active_ids = {
            record.track_id
            for record in table
            if record.start_frame <= frame <= record.end_frame
        }
        observed_ids = {int(region.label) for region in regionprops(mask)}
        if observed_ids != active_ids:
            details = []
            if active_ids - observed_ids:
                details.append(
                    "missing labels "
                    + ", ".join(str(v) for v in sorted(active_ids - observed_ids))
                )
            if observed_ids - active_ids:
                details.append(
                    "unexpected labels "
                    + ", ".join(str(v) for v in sorted(observed_ids - active_ids))
                )
            raise ValueError(
                f"CTC table and mask {path.name} disagree: " + "; ".join(details)
            )
        frames.append(mask)

    stack = np.stack(frames, axis=0)
    man_track = pd.DataFrame(
        [
            [
                record.track_id,
                record.start_frame,
                record.end_frame,
                record.parent_track_id,
            ]
            for record in table
        ]
    )
    tracks, graph = ctc_to_napari_tracks(stack, man_track)

    header = ["track_id", "t", "y", "x"]
    if stack.ndim == 4:
        header = ["track_id", "t", "z", "y", "x"]

    rows = [
        [int(entry[0]), int(entry[1]), *[float(v) for v in entry[2:]]]
        for entry in tracks
    ]
    rows.sort(key=lambda row: (row[0], row[1]))
    graph = {int(child): [int(v) for v in parents] for child, parents in graph.items()}
    return header, rows, graph, stack


def write_tracks(path: Path, header: list[str], rows: list[list[float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_graph(path: Path, graph: dict[int, list[int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(graph, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_napari_stack(
    path: Path,
    stack: "np.ndarray",
    header: list[str],
    rows: list[list[float | int]],
    graph: dict[int, list[int]],
    table: list[TrackRecord],
) -> None:
    """Write the tracked masks as one multi-page TIFF with the tracks embedded."""
    axes = "TYX" if stack.ndim == 3 else "TZYX"
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path,
        stack,
        compression="deflate",
        photometric="minisblack",
        metadata={
            "axes": axes,
            "trackastra": {
                "napari_tracks_columns": header,
                "napari_tracks": rows,
                "napari_graph": {
                    str(child): parents for child, parents in graph.items()
                },
                "ctc_track_table": [
                    [
                        record.track_id,
                        record.start_frame,
                        record.end_frame,
                        record.parent_track_id,
                    ]
                    for record in table
                ],
            },
        },
    )

    with tifffile.TiffFile(path) as handle:
        written = handle.asarray()
        metadata = handle.shaped_metadata[0]
    if written.shape != stack.shape or written.dtype != stack.dtype:
        raise RuntimeError(
            f"Napari stack changed while writing {path.name}: "
            f"{stack.shape}/{stack.dtype} -> {written.shape}/{written.dtype}"
        )
    if not np.array_equal(written, stack):
        raise RuntimeError(f"Pixel labels changed while writing {path.name}.")
    if metadata.get("axes") != axes:
        raise RuntimeError(
            f"Axis metadata was not stored in {path.name}; expected {axes}."
        )
    if len(metadata["trackastra"]["napari_tracks"]) != len(rows):
        raise RuntimeError(
            f"Track metadata was not stored completely in {path.name}."
        )


def write_ctc_directory(
    directory: Path,
    masks: list[tuple[int, Path]],
) -> None:
    """Write the validated masks under canonical CTC names."""
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise ValueError(f"CTC output directory is not empty: {directory}")
    for frame, mask in masks:
        destination = directory / f"man_track{frame:04d}{GALAXY_TIFF_SUFFIX}"
        shutil.copyfile(mask, destination)


def main() -> None:
    args = parse_args()
    table = read_track_table(args.track_table)
    masks = find_masks(args.masks_dir)
    header, rows, graph, stack = build_napari_outputs(table, masks)
    write_napari_stack(args.stack_out, stack, header, rows, graph, table)
    write_tracks(args.tracks_out, header, rows)
    write_graph(args.graph_out, graph)
    write_ctc_directory(args.ctc_dir, masks)


if __name__ == "__main__":
    main()
