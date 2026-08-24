#!/usr/bin/env python3
"""Validate and stage an ordered Galaxy collection for clinker."""

import argparse
import os
import re
from pathlib import Path

from Bio import SeqIO


def safe_name(label, fallback):
    """Return a readable, filesystem-safe name for a collection element."""
    stem = re.sub(
        r"\.(?:gbk|gb|genbank|gbf|gbff)$",
        "",
        label.strip(),
        flags=re.IGNORECASE,
    )
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem)
    stem = stem.strip("._-")
    return stem or fallback


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory in which ordered GBFF symlinks will be created.",
    )
    parser.add_argument(
        "--genbank",
        action="append",
        nargs=2,
        metavar=("PATH", "LABEL"),
        required=True,
        help="GenBank path and Galaxy collection element identifier.",
    )
    return parser.parse_args()


def validate_genbank(source, label):
    """Validate one input and return its record and CDS counts."""
    source = Path(source)
    if not source.is_file():
        raise SystemExit(f"Input for collection element {label!r} does not exist: {source}")

    records = 0
    cds_features = 0
    try:
        for record in SeqIO.parse(source, "genbank"):
            records += 1
            cds_features += sum(feature.type == "CDS" for feature in record.features)
    except (OSError, UnicodeError, ValueError) as error:
        raise SystemExit(
            f"Collection element {label!r} is not valid GenBank data: {error}"
        ) from error

    if records == 0:
        raise SystemExit(
            f"Collection element {label!r} contains no GenBank records"
        )
    if cds_features == 0:
        raise SystemExit(
            f"Collection element {label!r} contains no CDS features"
        )
    return records, cds_features


def stage_inputs(inputs, input_dir):
    """Validate inputs and stage them in collection order."""
    if len(inputs) < 2:
        raise SystemExit("clinker requires at least two GenBank collection elements")

    input_dir = Path(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)
    used_names = set()

    for index, (source, label) in enumerate(inputs, start=1):
        validate_genbank(source, label)

        base_name = safe_name(label, f"plasmid_{index}")
        staged_name = base_name
        suffix = 2
        while staged_name in used_names:
            staged_name = f"{base_name}_{suffix}"
            suffix += 1
        used_names.add(staged_name)

        destination_dir = input_dir / f"{index:06d}"
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination = destination_dir / f"{staged_name}.gbff"
        os.symlink(Path(source).resolve(), destination)


def main():
    args = parse_args()
    stage_inputs(args.genbank, args.input_dir)


if __name__ == "__main__":
    main()
