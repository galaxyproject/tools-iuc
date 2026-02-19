import argparse
import bz2
import gzip

from Bio import SeqIO


def get_opener(gbk_filename):
    """Determines the appropriate opener for a given file, supporting
    bzip2, gzip, or standard open.
    """
    try:
        bz2.open(gbk_filename).read(1)
        return bz2.open
    except OSError:
        pass
    try:
        gzip.open(gbk_filename).read(1)
        return gzip.open
    except OSError:
        return open


def main():
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to FASTA format. "
                    "Supports gzip and bzip2 compressed files."
    )
    parser.add_argument(
        "genbank_file",
        help="GenBank input file. Can be compressed with gzip or bzip2"
    )
    parser.add_argument(
        "fasta_file",
        help="FASTA output dataset"
    )
    parser.add_argument(
        "--remove_version", action="store_true",
        help="Remove version number from NCBI formatted accession numbers. "
             "For example, this converts 'B000657.2' to 'B000657'."
    )
    args = parser.parse_args()

    gbk_open = get_opener(args.genbank_file)
    with gbk_open(args.genbank_file, 'rt') as input_handle, \
            open(args.fasta_file, 'w') as output_handle:
        for seq_record in SeqIO.parse(input_handle, 'genbank'):
            if args.remove_version:
                seq_id = seq_record.id.split('.')[0]
            else:
                seq_id = seq_record.id
            print(f'Writing FASTA record: {seq_id}')
            output_handle.write(f'>{seq_id}\n')
            output_handle.write(f'{seq_record.seq}\n')


if __name__ == "__main__":
    main()
