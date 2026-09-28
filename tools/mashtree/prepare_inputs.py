#!/usr/bin/env python3
"""Validate and stage a Galaxy FASTA collection for Mashtree."""

import argparse
import gzip
import os
import re
from pathlib import Path


# Mashtree recursively strips these suffixes when deriving sample names.
MASHTREE_STRIPPED_EXTENSIONS = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
    ".fasta",
    ".fna",
    ".faa",
    ".mfa",
    ".fas",
    ".fsa",
    ".fa",
    ".gbank",
    ".genbank",
    ".gbk",
    ".gbs",
    ".gbf",
    ".gb",
    ".embl",
    ".ebl",
    ".emb",
    ".dat",
    ".swiss",
    ".sp",
    ".msh",
)


def safe_name(label, fallback):
    """Return a readable Mashtree-safe stem from a collection identifier."""
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", label.strip())
    stem = stem.strip("._-") or fallback

    lower_stem = stem.lower()
    for suffix in MASHTREE_STRIPPED_EXTENSIONS:
        if lower_stem.endswith(suffix):
            stem = stem[:-len(suffix)] + "_" + stem[-len(suffix) + 1:]
            break

    return stem


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory in which staged FASTA inputs will be created.",
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to write the list of staged input paths.",
    )
    parser.add_argument(
        "--fasta",
        action="append",
        nargs=2,
        metavar=("PATH", "LABEL"),
        required=True,
        help="FASTA path and Galaxy collection element identifier.",
    )
    return parser.parse_args()


def is_gzip(path):
    """Detect gzip content independently of the Galaxy dataset extension."""
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def validate_fasta(path, label, compressed):
    """Require at least one FASTA record containing sequence data."""
    opener = gzip.open if compressed else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            saw_header = False
            saw_sequence = False
            in_record = False
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if len(line) == 1:
                        raise SystemExit(
                            f"Collection element {label!r} contains an empty "
                            "FASTA header"
                        )
                    saw_header = True
                    in_record = True
                    continue
                if not in_record:
                    raise SystemExit(
                        f"Collection element {label!r} is not valid FASTA: "
                        "sequence data occurs before the first FASTA header"
                    )
                saw_sequence = True
    except (OSError, UnicodeError) as error:
        raise SystemExit(
            f"Could not read FASTA collection element {label!r}: {error}"
        ) from error

    if not saw_header:
        raise SystemExit(f"Collection element {label!r} contains no FASTA records")
    if not saw_sequence:
        raise SystemExit(f"Collection element {label!r} contains no sequence data")


def stage_inputs(inputs, input_dir, manifest_path):
    """Validate inputs and stage them without changing sequence content."""
    if len(inputs) < 2:
        raise SystemExit("Mashtree requires at least two FASTA collection elements")

    input_dir = Path(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)

    seen_labels = set()
    staged_names = {}
    manifest_paths = []

    for index, (source, label) in enumerate(inputs, start=1):
        source = Path(source)
        if not source.is_file():
            raise SystemExit(
                f"Input for collection element {label!r} does not exist: {source}"
            )
        if label in seen_labels:
            raise SystemExit(
                f"Collection element identifiers must be unique: {label!r}"
            )
        seen_labels.add(label)

        compressed = is_gzip(source)
        validate_fasta(source, label, compressed)

        stem = safe_name(label, f"sample_{index}")
        if stem in staged_names:
            other = staged_names[stem]
            raise SystemExit(
                "Collection element identifiers become ambiguous after filename "
                f"normalization: {other!r} and {label!r} both become {stem!r}"
            )
        staged_names[stem] = label

        suffix = ".fasta.gz" if compressed else ".fasta"
        destination_dir = input_dir / f"{index:06d}"
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination = destination_dir / f"{stem}{suffix}"
        os.symlink(source.resolve(), destination)
        manifest_paths.append(destination.absolute())

    with open(manifest_path, "w", encoding="utf-8") as handle:
        for path in manifest_paths:
            handle.write(f"{path}\n")


def main():
    args = parse_args()
    stage_inputs(args.fasta, Path(args.input_dir), Path(args.manifest))


if __name__ == "__main__":
    main()
