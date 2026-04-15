#!/usr/bin/env python3

#########################################################################################
# This script detect circular contigs by looking for exact identical k-mer at the two
# ends of the sequences provided in fasta file. In order to be able to predict genes
# spanning the origin of circular contigs, the first 1,000 nucleotides of each circular
# contigs are duplicated and added at the contig's end.
#
# Inspired by Simon Roux work for Metavir2 (2014) and Corentin Hochart work in PlasSuite
#
#########################################################################################

import argparse
import logging
from pathlib import Path

from Bio import SeqIO

log_levels = {
    0: logging.CRITICAL,
    1: logging.ERROR,
    2: logging.WARN,
    3: logging.INFO,
    4: logging.DEBUG,
}
logging.basicConfig(level=log_levels[3])
logger = logging.getLogger()


def setup_logger(verbosity: int) -> None:
    """
    Configure the logger based on verbosity level.

    :param verbosity: verbosity level
    """
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=log_levels.get(verbosity, logging.INFO),
    )


def find_occurrences(s, substring) -> list:
    """
    Find all starting positions of a substring in a string

    :param s: String to be searched
    :param substring: Substring to search in s
    """
    return [i for i in range(len(s)) if s.startswith(substring, i)]


def is_circular(sequence, length, pos) -> bool:
    """
    Determines if a sequence is circular by comparing segments starting at `start_pos`.

    A sequence is considered circular if the `length` elements at the beginning of the sequence match the `length`
    elements starting at `start_pos` in the sequence. This is useful for detecting repeating patterns or cycles
    in sequences.

    :param sequence: The input sequence
    :param length: The number of elements to compare for circularity.
    :param pos: The starting index in the sequence to begin the comparison.

    :return bool: True if circular, False otherwise
    """
    for i in range(length):
        if sequence[i] != sequence[pos + i]:
            return False
    return True


def check_circularity(seq_record, subseq_length=10) -> int:
    """
    Process a single sequence to detect circularity and return the overlap length if circular.

    :param seq_record: SeqRecord object
    :param subseq_length: Length of 3' fragment to check on the 5' end

    :return: overlap length if circular, 0 otherwise
    """
    seq_len = len(seq_record)

    if seq_len < subseq_length:
        logging.error(f"Sequence too short ({seq_len}bp): {seq_record.id}")
        return 0

    begin = "".join(seq_record[:subseq_length])
    end = "".join(seq_record[subseq_length:])
    positions = [x + subseq_length for x in find_occurrences(end, begin)]

    for pos in positions:
        overlap_length = seq_len - pos
        if is_circular(seq_record, overlap_length, pos):
            return overlap_length
    return 0


def extend_sequence(
    seq_record,
    overlap_length,
    duplication_length=1000,
):
    """
    Extends the 5' end of a sequence by duplicating a fragment from the 3' end.

    This function is useful for simulating circular sequences by extending the 5' end
    with a fragment from the 3' end, based on the specified `overlap_length` and `duplication_length`.

    :param seq_record: The input sequence record to be extended.
    :param overlap_length: The length of the overlapping segment that was previously identified as circular.
    :param duplication_length: The length of the 3' end fragment to duplicate and add to the 5' end.

    :return: The modified sequence record with the extended 5' end.
    """
    # Remove the overlapping segment from the 3' end
    modified_seq = seq_record.seq[: len(seq_record.seq) - overlap_length]
    # Duplicate the first `duplication_length` nucleotides from the original sequence
    # and append them to the 5' end of the modified sequence
    if len(modified_seq) < duplication_length:
        # If the modified sequence is shorter than `duplication_length`,
        # duplicate the entire modified sequence
        extension = modified_seq
    else:
        # Otherwise, duplicate the first `duplication_length` nucleotides
        extension = seq_record.seq[:duplication_length]
    # Combine the modified sequence with the duplicated fragment
    extended_seq = modified_seq + extension
    # Update the sequence in the SeqRecord object
    seq_record.seq = extended_seq
    return seq_record


def detect_circular(
    fasta_in,
    fasta_out,
    id_out,
    subseq_length=10,
    duplication_length=1000,
):
    """
    Detect and process circular sequences in a FASTA file.

    This function reads sequences from `fasta_in`, checks for circularity,
    extends circular sequences, and writes the results to `fasta_out` and `id_out`.

    :param fasta_in: Path to the input FASTA file.
    :param fasta_out: Path to the output FASTA file for extended circular sequences.
    :param id_out: Path to the output file for recording IDs of circular sequences.
    :param subseq_length: Length of the 3' fragment to check for circularity.
    :param duplication_length: Length of the 3' fragment to duplicate and add to the 5' end.
    """
    records = []
    ids = []
    try:
        with fasta_in.open("r") as fasta_in_f:
            for seq_record in SeqIO.parse(fasta_in_f, "fasta"):
                overlap_length = check_circularity(
                    seq_record,
                    subseq_length=subseq_length,
                )
                if overlap_length > 0:
                    records.append(
                        extend_sequence(
                            seq_record,
                            overlap_length,
                            duplication_length,
                        )
                    )
                    ids.append(seq_record.id)
    except Exception as e:
        logging.error(f"Error processing {fasta_in}: {e}")
        raise

    if not records:
        logging.warning("Warning: No circular sequences found.")

    try:
        with fasta_out.open("w") as fasta_out_f:
            SeqIO.write(records, fasta_out_f, "fasta")
        with id_out.open("w") as id_out_f:
            id_out_f.write("\n".join(ids) + "\n")
    except IOError as e:
        logging.error(f"Error writing output files: {e}")
        raise


def main():
    """
    Main function to detect circular contigs in a FASTA file.

    This function parses command-line arguments, launches function to read the input
    FASTA file, process each sequence to detect circular contigs, and generate the
    output files.
    """
    parser = argparse.ArgumentParser(
        description="""
            Detect circular contigs by looking for exact identical subsequences at the two
            ends of the sequences provided in a FASTA file and output the circular contigs
            extended on 5' end by duplication of the first nucleotides on 3' end to be able
            to predict genes spanning the origin of circular contigs.
        """
    )
    parser.add_argument("--fasta-in", required=True, help="Input FASTA file")
    parser.add_argument(
        "--subseq-length",
        type=int,
        default=10,
        help="Length of 3' fragment to check on the 5' end (default: 10)",
    )
    parser.add_argument(
        "--duplication-length",
        type=int,
        default=1000,
        help="Length of the 3' end fragment to duplicate and add on the 5' end (default: 1000)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=3,
        choices=log_levels.keys(),
        help="Verbosity level (0=CRITICAL, 1=ERROR, 2=WARN, 3=INFO, 4=DEBUG)",
    )
    parser.add_argument(
        "--fasta-out",
        required=True,
        help="Output FASTA file with extended circular contigs",
    )
    parser.add_argument(
        "--id-out", required=True, help="Output TXT file with circular sequence IDs"
    )

    args = parser.parse_args()
    setup_logger(args.verbose)

    logging.info("Starting script execution.")
    detect_circular(
        Path(args.fasta_in),
        Path(args.fasta_out),
        Path(args.id_out),
        subseq_length=args.subseq_length,
        duplication_length=args.duplication_length,
    )
    logging.info("Script execution completed.")


if __name__ == "__main__":
    main()
