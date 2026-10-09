#!/usr/bin/env python3
"""
Remove terminal stop codons from coding sequences.

Trim all terminal stop codons from sequences in a FASTA file, using a chosen
NCBI genetic code (translation table). Leave non-stop terminal codons alone.
If any INTERNAL, in-frame stop codon is found, exit with an error.

This tool is designed as a preprocessing step for tools like cawlign and HyPhy
that do not permit internal stop codons in their input sequences.

Requires: Biopython
"""

import argparse
import sys

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq


def load_table(table_arg):
    """Return a DNA codon table from an NCBI table id (int) or name (str)."""
    if table_arg is None:
        return CodonTable.unambiguous_dna_by_id[1]  # Standard
    # try as integer id
    try:
        tid = int(table_arg)
        return CodonTable.unambiguous_dna_by_id[tid]
    except (ValueError, KeyError):
        pass
    # try as name
    try:
        return CodonTable.unambiguous_dna_by_name[table_arg]
    except KeyError:
        # Build a helpful hint list
        valid_ids = sorted(CodonTable.unambiguous_dna_by_id.keys())
        valid_names = sorted(CodonTable.unambiguous_dna_by_name.keys())
        sys.stderr.write(
            f"ERROR: Unknown genetic code '{table_arg}'.\n"
            f"Try an NCBI table id (e.g., 1) or one of these names:\n"
            f"  {', '.join(valid_names)}\n"
            f"(Valid ids include: {', '.join(map(str, valid_ids))})\n"
        )
        sys.exit(2)


def trim_terminal_stops_and_validate(record, stop_codons, check_internal=True):
    """
    Remove ALL trailing stop codons (0+ at the end).

    If check_internal is True and any internal in-frame stop codon exists
    (excluding the trailing block), exit with an error message.

    Ignore a terminal codon that is not a stop codon.
    """
    # Work with DNA letters; treat any RNA U as T
    seq_str = str(record.seq).upper().replace("U", "T")

    # Count how many full codons sit at the end that are stops
    idx = len(seq_str)
    trailing_stops = 0
    while idx >= 3:
        codon = seq_str[idx - 3:idx]
        if codon in stop_codons:
            trailing_stops += 1
            idx -= 3
        else:
            break

    # Scan for INTERNAL stops: all complete codons up to (but not including)
    # the trailing stop block (and ignoring any trailing partial codon).
    if check_internal:
        scan_end = (idx // 3) * 3  # only complete codons
        for pos in range(0, scan_end, 3):
            codon = seq_str[pos:pos + 3]
            if codon in stop_codons:
                sys.stderr.write(
                    f"ERROR: Found an internal stop codon in sequence "
                    f"'{record.id}' at position {pos}.\n"
                    f"Tools like HyPhy and cawlign do not permit internal "
                    f"stop codons. Please review your input sequences.\n"
                )
                sys.exit(2)

    # Finally, remove the trailing stop codons (if any)
    if trailing_stops > 0:
        seq_str = seq_str[:idx]

    # Leave sequences with non-stop terminal codons unchanged by design
    return Seq(seq_str)


def main():
    ap = argparse.ArgumentParser(
        description="Remove all terminal stop codons from a FASTA, using a "
                    "chosen genetic code. Optionally fail if any internal "
                    "in-frame stop codon is present."
    )
    ap.add_argument("-i", "--input", required=True, help="Input FASTA file")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA file")
    ap.add_argument(
        "-t", "--table",
        help="NCBI translation table id (e.g., 1) or name "
             "(e.g., 'Vertebrate Mitochondrial'). Default: 1 (Standard)."
    )
    ap.add_argument(
        "--no-check-internal",
        action="store_true",
        help="Do not check for internal stop codons (only remove terminal)."
    )
    args = ap.parse_args()

    table = load_table(args.table)
    stop_codons = set(table.stop_codons)  # e.g., {'TAA','TAG','TGA'} for Standard

    check_internal = not args.no_check_internal

    records_out = []
    for rec in SeqIO.parse(args.input, "fasta"):
        new_seq = trim_terminal_stops_and_validate(rec, stop_codons, check_internal)
        rec.seq = new_seq
        records_out.append(rec)

    SeqIO.write(records_out, args.output, "fasta")


if __name__ == "__main__":
    main()
