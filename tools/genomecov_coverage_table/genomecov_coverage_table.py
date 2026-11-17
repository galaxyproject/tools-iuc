#!/usr/bin/env python3
"""Tabulate percent genome coverage from bedtools genomecov histograms.

The config is a JSON list of columns, in output order::

    [{"label": "dustmasker",
      "histograms": [{"identifier": "sampleA", "path": "/data/a.tabular"}, ...]},
     ...]

Each histogram is ``bedtools genomecov`` output in the default histogram mode. The
script sums column 5 (fraction) of the ``genome`` rows with depth >= ``--min-depth``
and reports it as a percent. A histogram with no such row reports 0.

The output header is ``Sample`` and the labels. There is one row per identifier, in
the order of the first column. Every column must have the same identifiers.
"""

import argparse
import json
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path

__version__ = "1.0.0"

LOG = logging.getLogger("genomecov_coverage_table")

FORBIDDEN_LABEL_CHARS = ("\t", "\n", "\r")


class InputError(Exception):
    """A problem with the user's inputs, reported without a traceback."""


@dataclass(frozen=True)
class Histogram:
    identifier: str
    path: Path


@dataclass(frozen=True)
class Column:
    label: str
    histograms: list[Histogram]


def require_str(value: object, *, what: str) -> str:
    if not isinstance(value, str):
        raise InputError(f"{what} must be a string, got {type(value).__name__}")
    return value


def parse_histogram(entry: object, *, label: str) -> Histogram:
    if not isinstance(entry, dict):
        raise InputError(f"column '{label}': each histogram must be an object")
    if "identifier" not in entry or "path" not in entry:
        raise InputError(
            f"column '{label}': each histogram needs 'identifier' and 'path'"
        )
    return Histogram(
        identifier=require_str(
            entry["identifier"], what=f"column '{label}' identifier"
        ),
        path=Path(require_str(entry["path"], what=f"column '{label}' path")),
    )


def parse_column(entry: object) -> Column:
    if not isinstance(entry, dict):
        raise InputError("each column must be an object")
    if "label" not in entry or "histograms" not in entry:
        raise InputError("each column needs 'label' and 'histograms'")
    label = require_str(entry["label"], what="label")
    histograms = entry["histograms"]
    if not isinstance(histograms, list):
        raise InputError(f"column '{label}': 'histograms' must be a list")
    return Column(
        label=label, histograms=[parse_histogram(h, label=label) for h in histograms]
    )


def load_columns(config_path: Path) -> list[Column]:
    try:
        data = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise InputError(f"config {config_path} is not valid JSON: {e}") from e
    if not isinstance(data, list) or not data:
        raise InputError("config must be a non-empty list of columns")
    return [parse_column(entry) for entry in data]


def validate_labels(columns: list[Column]) -> None:
    seen: set[str] = set()
    for column in columns:
        if column.label.strip() == "":
            raise InputError("a column label is empty")
        if any(char in column.label for char in FORBIDDEN_LABEL_CHARS):
            raise InputError(f"column label {column.label!r} contains a tab or newline")
        if column.label in seen:
            raise InputError(f"duplicate column label '{column.label}'")
        seen.add(column.label)


def sample_order(columns: list[Column]) -> list[str]:
    """Return the first column's identifiers after checking every column has the same set."""
    for column in columns:
        identifiers = [h.identifier for h in column.histograms]
        if len(set(identifiers)) != len(identifiers):
            raise InputError(f"column '{column.label}' has duplicate identifiers")
    first = columns[0]
    order = [h.identifier for h in first.histograms]
    expected = set(order)
    for column in columns[1:]:
        found = {h.identifier for h in column.histograms}
        if found == expected:
            continue
        missing = ", ".join(sorted(expected - found)) or "none"
        extra = ", ".join(sorted(found - expected)) or "none"
        raise InputError(
            f"column '{column.label}' identifiers differ from column '{first.label}': "
            f"missing {missing}; extra {extra}"
        )
    return order


def covered_percent(path: Path, *, min_depth: int) -> float:
    """Return 100 x the summed ``genome`` fraction at depth >= min_depth."""
    fractions: list[float] = []
    genome_rows = 0
    with path.open(encoding="utf-8") as fh:
        for line_number, line in enumerate(fh, start=1):
            cols = line.rstrip("\n").split("\t")
            if cols[0] != "genome":
                continue
            if len(cols) < 5:
                raise InputError(
                    f"{path} line {line_number}: genome row has fewer than 5 columns"
                )
            try:
                depth = int(cols[1])
                fraction = float(cols[4])
            except ValueError as e:
                raise InputError(
                    f"{path} line {line_number}: bad genome row: {e}"
                ) from e
            genome_rows += 1
            if depth >= min_depth:
                fractions.append(fraction)
    if genome_rows == 0:
        raise InputError(
            f"{path} has no 'genome' rows; is it bedtools genomecov histogram output?"
        )
    return 100.0 * math.fsum(fractions)


def format_percent(value: float, *, decimals: int) -> str:
    return f"{round(value, decimals):.{decimals}f}"


def write_table(
    columns: list[Column], *, out: Path, min_depth: int, decimals: int
) -> None:
    samples = sample_order(columns)
    paths = [{h.identifier: h.path for h in column.histograms} for column in columns]
    lines = ["\t".join(["Sample", *(column.label for column in columns)])]
    for sample in samples:
        values = [
            format_percent(
                covered_percent(by_id[sample], min_depth=min_depth), decimals=decimals
            )
            for by_id in paths
        ]
        lines.append("\t".join([sample, *values]))
    out.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
    LOG.info("wrote %d samples x %d columns to %s", len(samples), len(columns), out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--config", required=True, type=Path, help="JSON column config")
    parser.add_argument(
        "--min-depth", type=int, default=1, help="minimum depth counted as covered"
    )
    parser.add_argument(
        "--decimals", type=int, default=2, help="decimal places in the output"
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="output tabular file"
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main() -> None:
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr
    )
    args = build_parser().parse_args()
    if args.min_depth < 1:
        LOG.error("--min-depth must be at least 1")
        raise SystemExit(1)
    if args.decimals < 0:
        LOG.error("--decimals must not be negative")
        raise SystemExit(1)
    try:
        columns = load_columns(args.config)
        validate_labels(columns)
        write_table(
            columns, out=args.output, min_depth=args.min_depth, decimals=args.decimals
        )
    except InputError as e:
        LOG.error("%s", e)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
