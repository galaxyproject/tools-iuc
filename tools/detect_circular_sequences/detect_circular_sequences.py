#!/usr/bin/env python3

#########################################################################################
# This script detect circular contigs by looking for exact identical k-mer at the two
# ends on a cadre sequence of the sequences prodvide in fasta file. In order to be able
# to predict genes spanning the orgin of circular contigs, the first 1,000 nucleotides
# of each circular contigs are dulicated and added at the contig's end.
#
# Inspired by Simon Roux work for Metavir2 (2014) and Corentin Hochart work in PlasSuite
#
#########################################################################################

import argparse
import re
import sys
import tempfile
import textwrap
from pathlib import Path


def error(message):
    """
    Print an error message to stderr and exit the program with a non-zero status.

    Args:
        message (str): The error message to display.
    """
    print(f"{sys.argv[0]} (error): {message}. Execution halted.", file=sys.stderr)
    sys.exit(1)


def warning(message, verbose):
    """
    Print a warning message to stderr if verbose mode is enabled.

    Args:
        message (str): The warning message to display.
        verbose (bool): If True, the message will be printed.
    """
    if verbose:
        print(f"{sys.argv[0]} (info): {message}", file=sys.stderr)


def fasta_format(seq):
    """
    Format sequence into lines of 60 characters each.

    Args:
        seq (str): sequence to format.
    """
    return textwrap.wrap(seq, width=60, break_on_hyphens=False)


def one_line_fasta(input_fp, output_fp):
    """
    Convert FASTA file to a format with sequences on single lines.

    Args:
        input_fp (Path): path to input FASTA file
        output_fp (Path): path to output FASTA file
    """
    with input_fp.open("r") as infile, output_fp.open("w") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(f"\n{line}")  # Newline before header
            else:
                outfile.write(line.rstrip("\n"))  # Remove newline and concatenate
        outfile.write("\n")  # Final newline (like END in awk)


def find_kmer_occurrences(begin, end):
    """
    Find all starting positions of 'begin' in 'end'.

    Args:
        begin ():
        end ():
    """
    pattern = re.compile(re.escape(begin))
    return [match.start() + len(begin) for match in pattern.finditer(end)]


def is_circular(line_chars, scale, pos):
    """
    Check if the sequence is circular by comparing segments.

    Args:
        line_chars (list): Sequence characters
        scale ( ): Starting k-mer
        pos ():

    Returns:
        bool: True if circular, False otherwise
    """
    for i in range(scale):
        if line_chars[i] != line_chars[pos + i]:
            return False
    return True


def process_sequence(
    line,
    header,
    verbose,
    kmer_length=10,
    cadre_length=0,
    duplicate_nucleotides=1000,
):
    """
    Process a single sequence to detect circularity and modify if needed.

    Args:
        line ():
        header ():
        verbose ():
        kmer_length ():
        cadre_length ():
        duplicate_nucleotides ():

    Returns:
        : True if circular, False otherwise
    """
    try:
        line_chars = list(line)
        seq_len = len(line_chars)

        if seq_len < kmer_length:
            warning(f"Short sequence ({seq_len}bp): {header}", verbose)
            return None, None

        # Determine cadre length
        if cadre_length == 0 or cadre_length > seq_len - kmer_length:
            cadre_length = seq_len - kmer_length

        # Extract begin and end
        begin = "".join(line_chars[:kmer_length])
        end_part = line_chars[-cadre_length:]  # [::-1]
        end_str = "".join(end_part)

        # Find all positions where 'begin' appears in the reversed end part
        end_positions = find_kmer_occurrences(begin, end_str)

        if not end_positions:
            return None, None

        # Check for circularity at each position
        status = False
        for pos in end_positions:
            scale = len(line_chars) - pos
            if is_circular(line_chars, scale, pos):
                status = True

        if not status:
            return None, None

        # Modify the sequence
        modified_seq = line_chars[: len(line_chars) - scale]

        if len(modified_seq) < duplicate_nucleotides:
            modified_seq += modified_seq
        else:
            modified_seq += line_chars[:duplicate_nucleotides]

        return header, "".join(modified_seq)
    except Exception as e:
        error(f"Error processing sequence {header}: {e}")


def detect_circular(
    fasta_in,
    fasta_out,
    id_out,
    kmer_length=10,
    cadre_length=0,
    duplicate_nucleotides=1000,
    verbose=False,
):
    """
    Detect and process circular sequences.

    Args:
        fasta_in (Path): Input FASTA file
        fasta_out (Path): Output FASTA file with modifications
        id_out (Path): File to record identifiers of circular sequences
        kmer_length (int): Length of k-mer to search for
        cadre_length (int): Length of sequence end to inspect
        duplicate_nucleotides (int): Number of nucleotides to duplicate at the end
        verbose (bool): Enable verbose output
    """
    tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    one_line_fasta(fasta_in, Path(tmp_file.name))

    with Path(tmp_file.name).open("r") as infile, fasta_out.open(
        "w"
    ) as fasta_output, id_out.open("w") as id_output:

        sequence = ""
        header = ""
        for line in infile.readlines():
            # print(line)
            # Read header
            if line.startswith(">"):
                header = line[1:].strip()
                sequence = ""
                continue
            else:
                sequence = line.strip()

            # Process the sequence
            header, modified_seq = process_sequence(
                sequence,
                header,
                verbose=verbose,
                kmer_length=kmer_length,
                cadre_length=cadre_length,
                duplicate_nucleotides=duplicate_nucleotides,
            )

            if modified_seq:
                # Write to output files
                id_output.write(f"{header}\n")

                fasta_output.write(f">{header}\n")
                formatted = fasta_format(modified_seq)
                for line in formatted:
                    fasta_output.write(f"{line}\n")

    # Clean up temporary file
    Path(tmp_file.name).unlink()


def main():
    """
    Main function to detect circular contigs in a FASTA file.

    This function parses command-line arguments, reads the input FASTA file,
    processes each sequence to detect circular contigs, and prints the results.
    It handles both verbose and non-verbose modes, and can optionally print only circular sequences.

    Command-line arguments:
        --fasta-in: Path to the input FASTA file (required).
        --kmer-length: Length of the k-mer used to identify circular sequences (default: 10).
        --cadre: Length of the fragment at the sequence 5' end to inspect for k-mer identity.
                 If 0, the entire sequence is screened (default: 0).
        --only_circular: If set, only circular sequences are printed.
        --verbose: If set, warning messages are printed during execution.
        --output: Path to output file

    The function processes each sequence in the FASTA file, checks for circularity,
    and prints the results in FASTA format, with circular sequences marked in the header.
    """
    parser = argparse.ArgumentParser(
        description="Detect circular contigs by k-mer matching."
    )
    parser.add_argument("--fasta-in", required=True, help="Input FASTA file")
    parser.add_argument(
        "--kmer-length", type=int, default=10, help="Length of k-mer (default: 10)"
    )
    parser.add_argument(
        "--cadre-length",
        type=int,
        default=0,
        help="Inspect fragment length (default: 0)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--fasta-out", required=True, help="Output FASTA file")
    parser.add_argument(
        "--id-out", required=True, help="File to write circular sequence IDs"
    )

    args = parser.parse_args()

    warning("Starting script execution.", args.verbose)
    detect_circular(
        Path(args.fasta_in),
        Path(args.fasta_out),
        Path(args.id_out),
        kmer_length=args.kmer_length,
        cadre_length=args.cadre_length,
        verbose=args.verbose,
    )
    warning("Script execution completed.", args.verbose)


if __name__ == "__main__":
    main()
