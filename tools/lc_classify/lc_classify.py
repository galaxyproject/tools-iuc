#!/usr/bin/env python
"""Classify low-complexity intervals of a FASTA by their repeat unit.

Reads a FASTA (plain or gzip) and a BED of intervals on it (only the first
three columns are used) and writes BED6 to stdout, one or more rows per
interval, in input BED order:

  name   = repeat-unit signature: 'polyA', 'polyT', ... (period 1),
           '(AT)n', '(AAT)n', ... (period 2-6), or 'lc' (anything else)
  score  = purity x 1000, the fraction of bases that match the best periodic
           consensus, 0-1000
  strand = '.'

Sequence is uppercased before classification. Each interval is split at runs
of N, and each non-N piece is classified on its own. An interval that is all N
gives no output row.

The FASTA is streamed one record at a time so that memory stays at one
sequence plus the BED, not the whole assembly. The BED may be in any order.

A BED line that names a sequence not in the FASTA, or an interval outside its
sequence, is a fatal error. The script then writes nothing and exits 1. The
earliest offending BED line is reported.
"""

import argparse
import gzip
import logging
import re
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TextIO

__version__ = "1.0.0"

MAX_PERIOD = 6
# Purity at which the smallest period is accepted without trying longer ones.
GOOD_PURITY = 0.85
# Purity below which the interval is called 'lc'.
MIN_PURITY = 0.60
GZIP_MAGIC = b"\x1f\x8b"


class InputMismatchError(Exception):
    """The BED and the FASTA do not fit together."""


@dataclass(frozen=True)
class Interval:
    """One BED interval and the output slot that its rows go to."""

    line: int
    start: int
    end: int
    slot: int


@cache
def non_gap_pattern() -> re.Pattern[str]:
    """Match maximal runs of non-gap bases."""
    return re.compile(r"[^N]+")


def is_gzip(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(2) == GZIP_MAGIC


def open_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed file for text reading."""
    if is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (name, sequence) one record at a time."""
    name: str | None = None
    buf: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def classify(seq: str) -> tuple[str, int]:
    """Return (signature, score) for an uppercase sequence without N."""
    length = len(seq)
    if length == 0:
        return "lc", 0
    best_frac = 0.0
    best_unit = ""
    for period in range(1, min(MAX_PERIOD, length) + 1):
        unit, match = periodic_consensus(seq, period)
        frac = match / length
        if frac > best_frac + 1e-9:
            best_frac, best_unit = frac, unit
        if frac >= GOOD_PURITY:
            # The smallest period that explains the interval well wins.
            best_frac, best_unit = frac, unit
            break
    score = min(1000, max(0, round(best_frac * 1000)))
    if best_frac < MIN_PURITY:
        return "lc", score
    if len(best_unit) == 1:
        return "poly" + best_unit, score
    return "(" + best_unit + ")n", score


def periodic_consensus(seq: str, period: int) -> tuple[str, int]:
    """Return the consensus unit of seq at period and the bases that match it."""
    cols: list[dict[str, int]] = [{} for _ in range(period)]
    for i, base in enumerate(seq):
        col = cols[i % period]
        col[base] = col.get(base, 0) + 1
    unit = ""
    match = 0
    for col in cols:
        # max() keeps the first base seen on a tie, as dicts keep insertion order.
        base = max(col, key=col.__getitem__)
        unit += base
        match += col[base]
    return unit, match


def parse_coordinate(field: str, *, line: int) -> int:
    if not field.isascii() or not field.isdigit():
        raise InputMismatchError(
            f"lc_classify: BED line {line}: start and end must be non-negative integers."
        )
    return int(field)


def read_bed(lines: Iterable[str]) -> tuple[dict[str, list[Interval]], int]:
    """Group BED intervals by sequence name.

    Returns the groups and the number of intervals. Header, track and browser
    lines, and lines with fewer than three columns, are skipped.
    """
    by_seq: dict[str, list[Interval]] = {}
    count = 0
    for line_no, line in enumerate(lines, 1):
        if line.startswith(("#", "track", "browser")):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 3:
            continue
        interval = Interval(
            line=line_no,
            start=parse_coordinate(fields[1], line=line_no),
            end=parse_coordinate(fields[2], line=line_no),
            slot=count,
        )
        by_seq.setdefault(fields[0], []).append(interval)
        count += 1
    return by_seq, count


def classify_interval(name: str, seq: str, interval: Interval) -> str:
    """Return the BED6 rows for one interval, split at runs of N."""
    start, end = interval.start, interval.end
    seg = seq[start:end].upper()
    rows = []
    for m in non_gap_pattern().finditer(seg):
        signature, score = classify(m.group())
        rows.append(
            f"{name}\t{start + m.start()}\t{start + m.end()}\t{signature}\t{score}\t.\n"
        )
    return "".join(rows)


def classify_intervals(
    *, fasta: Path, bed: Path, by_seq: dict[str, list[Interval]], count: int
) -> list[str]:
    """Return the BED6 text for each interval, in BED order.

    Raises InputMismatchError for the earliest BED line that cannot be
    classified. by_seq is consumed.
    """
    rows = [""] * count
    faults: list[tuple[int, str]] = []
    for name, seq in iter_fasta(fasta):
        if name not in by_seq:
            continue
        for interval in by_seq.pop(name):
            if interval.start >= interval.end or interval.end > len(seq):
                message = (
                    f"lc_classify: {bed} line {interval.line} has interval "
                    f"{interval.start}-{interval.end}, outside {name!r} (length {len(seq)})."
                )
                faults.append((interval.line, message))
                continue
            rows[interval.slot] = classify_interval(name, seq, interval)
    for seq_name, intervals in by_seq.items():
        line = intervals[0].line
        message = (
            f"lc_classify: {bed} line {line} names sequence {seq_name!r}, which is not in {fasta}. "
            "The BED and the FASTA do not describe the same assembly."
        )
        faults.append((line, message))
    if faults:
        raise InputMismatchError(min(faults)[1])
    return rows


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Classify low-complexity intervals of a FASTA by their repeat unit."
    )
    parser.add_argument("fasta", type=Path, help="FASTA file, plain or gzip")
    parser.add_argument(
        "bed",
        type=Path,
        help="BED file of intervals on the FASTA (BED3 or more columns)",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser.parse_args(argv)


def main(argv: list[str]) -> None:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    log = logging.getLogger("lc_classify")
    try:
        with args.bed.open(encoding="utf-8") as fh:
            by_seq, count = read_bed(fh)
        log.info("read %d intervals on %d sequences", count, len(by_seq))
        rows = classify_intervals(
            fasta=args.fasta, bed=args.bed, by_seq=by_seq, count=count
        )
    except InputMismatchError as e:
        raise SystemExit(str(e)) from None
    sys.stdout.writelines(rows)
    log.info("wrote %d rows", sum(row.count("\n") for row in rows))


if __name__ == "__main__":
    main(sys.argv[1:])
