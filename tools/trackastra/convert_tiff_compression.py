#!/usr/bin/env python3
"""Rewrite Trackastra CTC masks with Galaxy-readable lossless compression."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import tifffile

MASK_PATTERN = re.compile(r"^man_track[0-9]+\.tif$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Trackastra's ZSTD-compressed CTC TIFF masks to "
            "lossless Deflate TIFFs for Galaxy metadata and display support."
        )
    )
    parser.add_argument("source", type=Path, help="Trackastra CTC output directory")
    parser.add_argument(
        "destination",
        type=Path,
        help="Galaxy-compatible output directory",
    )
    return parser.parse_args()


def convert_file(source: Path, destination: Path) -> None:
    """Convert one TIFF and verify that shape, dtype, and pixels are unchanged."""
    array = tifffile.imread(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.tif")

    try:
        if array.ndim not in (2, 3):
            raise RuntimeError(
                f"Expected a 2D or 3D CTC mask in {source.name}; "
                f"found shape {array.shape}."
            )
        axes = "YX" if array.ndim == 2 else "ZYX"
        tifffile.imwrite(
            temporary,
            array,
            compression="deflate",
            photometric="minisblack",
            metadata={"axes": axes},
        )
        converted = tifffile.imread(temporary)
        if converted.shape != array.shape:
            raise RuntimeError(
                f"Shape changed while converting {source.name}: "
                f"{array.shape} -> {converted.shape}"
            )
        if converted.dtype != array.dtype:
            raise RuntimeError(
                f"Data type changed while converting {source.name}: "
                f"{array.dtype} -> {converted.dtype}"
            )
        if not np.array_equal(converted, array):
            raise RuntimeError(
                f"Pixel labels changed while converting {source.name}."
            )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    if not args.source.is_dir():
        raise ValueError(f"CTC output directory does not exist: {args.source}")
    if args.destination.exists() and any(args.destination.iterdir()):
        raise ValueError(
            f"Destination directory is not empty: {args.destination}"
        )

    masks = sorted(
        path
        for path in args.source.iterdir()
        if path.is_file() and MASK_PATTERN.fullmatch(path.name)
    )
    if not masks:
        raise ValueError(
            f"No Trackastra CTC mask TIFFs were found in {args.source}."
        )

    args.destination.mkdir(parents=True, exist_ok=True)
    for source in masks:
        convert_file(source, args.destination / source.name)


if __name__ == "__main__":
    main()
