#!/usr/bin/env python
"""Create Pling input files from multiple Galaxy FASTA datasets.

Pling takes a text manifest of FASTA paths rather than Galaxy's multi-dataset
input.
This helper sorts datasets by Galaxy element identifier, gives each one a stable
FASTA filename, writes the manifest, and optionally rewrites a topology TSV so
that its plasmid names match the generated FASTA stems.
"""

import argparse
import csv
import gzip
import os
import re
import shutil
from pathlib import Path

TOPOLOGY_FIELDS = ["plasmid", "topology"]
TOPOLOGY_VALUES = {"circular", "linear"}


def safe_name(label, fallback):
    """Return a filesystem-friendly FASTA stem derived from an input identifier."""
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", label.strip())
    stem = stem.strip("._-")
    return stem or fallback


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the Pling genomes_list manifest.",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory where prepared FASTA inputs are created.",
    )
    parser.add_argument(
        "--genome",
        nargs=2,
        action="append",
        metavar=("PATH", "LABEL"),
        required=True,
        help=(
            "Input FASTA path and Galaxy dataset identifier. "
            "May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--topology-input",
        help="Optional topology TSV with columns plasmid and topology.",
    )
    parser.add_argument(
        "--topology-output",
        help="Path to write the topology TSV remapped to generated names.",
    )
    args = parser.parse_args()

    if bool(args.topology_input) != bool(args.topology_output):
        raise SystemExit(
            "--topology-input and --topology-output must be supplied together"
        )
    return args


def is_gzip(path):
    """Detect gzip content independently of the Galaxy dataset extension."""
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def prepare_fasta(source, destination):
    """Materialize one uncompressed FASTA at the generated input path."""
    # Galaxy installations can preserve or transparently decompress compressed
    # datasets, so inspect the file content rather than trusting its extension.
    if is_gzip(source):
        try:
            with gzip.open(source, "rb") as compressed:
                with open(destination, "wb") as uncompressed:
                    shutil.copyfileobj(compressed, uncompressed)
        except (OSError, EOFError) as error:
            raise SystemExit(f"Could not decompress FASTA input: {source}") from error
    else:
        os.symlink(os.path.abspath(source), destination)


def prepare_inputs(genomes, input_dir):
    """Prepare sorted FASTAs and map dataset identifiers to generated stems."""
    input_dir.mkdir(parents=True, exist_ok=True)
    used_names = set()
    manifest_paths = []
    label_to_name = {}

    # Pling results can depend on genomes_list order. Sorting Galaxy dataset
    # identifiers here makes repeated runs with the same inputs reproducible.
    sorted_genomes = sorted(genomes, key=lambda genome: genome[1])
    for index, (source, label) in enumerate(sorted_genomes, start=1):
        if label in label_to_name:
            raise SystemExit(
                f"Input dataset identifiers must be unique: {label}"
            )
        base_name = safe_name(label, f"genome_{index}")
        name = base_name
        suffix = 2
        # Distinct labels can collapse to the same sanitized filename.
        while name in used_names:
            name = f"{base_name}_{suffix}"
            suffix += 1
        used_names.add(name)

        destination = input_dir / f"{name}.fasta"
        prepare_fasta(source, destination)
        manifest_paths.append(destination.absolute())

        label_to_name[label] = name

    return manifest_paths, label_to_name


def write_manifest(path, manifest_paths):
    """Write one prepared FASTA path per line for Pling."""
    with open(path, "w", encoding="utf-8") as handle:
        for manifest_path in manifest_paths:
            handle.write(f"{manifest_path}\n")


def read_topology(path, label_to_name):
    """Validate a topology TSV and remap labels to prepared FASTA stems."""
    with open(path, newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Topology TSV is empty")
        required = set(TOPOLOGY_FIELDS)
        missing = required.difference(reader.fieldnames)
        if missing:
            missing_names = ", ".join(sorted(missing))
            raise SystemExit(
                f"Topology TSV is missing required column(s): {missing_names}"
            )
        topology_by_label = {}
        for row in reader:
            plasmid = (row.get("plasmid") or "").strip()
            topology = (row.get("topology") or "").strip()
            if plasmid not in label_to_name:
                raise SystemExit(
                    "Topology TSV plasmid value does not match a selected "
                    f"dataset identifier: {plasmid}"
                )
            if plasmid in topology_by_label:
                raise SystemExit(
                    f"Topology TSV contains a duplicate plasmid value: {plasmid}"
                )
            if topology not in TOPOLOGY_VALUES:
                valid_values = ", ".join(sorted(TOPOLOGY_VALUES))
                raise SystemExit(
                    f"Invalid topology for {plasmid}: {topology}. "
                    f"Expected one of: {valid_values}"
                )
            topology_by_label[plasmid] = topology

    # Requiring every input avoids silently treating omitted linear plasmids as
    # circular, which is Pling's fallback when no topology entry is available.
    missing_plasmids = set(label_to_name).difference(topology_by_label)
    if missing_plasmids:
        missing_names = ", ".join(sorted(missing_plasmids))
        raise SystemExit(
            "Topology TSV is missing selected dataset(s): "
            f"{missing_names}"
        )

    # Pling matches topology IDs to the stems of the prepared FASTA paths.
    return [
        {
            "plasmid": generated_name,
            "topology": topology_by_label[label],
        }
        for label, generated_name in label_to_name.items()
    ]


def write_topology(path, rows):
    """Write the validated topology rows consumed by Pling."""
    with open(path, "w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(
            target,
            delimiter="\t",
            fieldnames=TOPOLOGY_FIELDS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    manifest_paths, label_to_name = prepare_inputs(
        args.genome,
        Path(args.input_dir),
    )
    write_manifest(args.output, manifest_paths)

    if args.topology_input:
        topology_rows = read_topology(args.topology_input, label_to_name)
        write_topology(args.topology_output, topology_rows)


if __name__ == "__main__":
    main()
