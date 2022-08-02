import argparse
import bz2
import gzip

from Bio import SeqIO


def get_opener(gbk_filename):
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


parser = argparse.ArgumentParser()
parser.add_argument(
    "genbank_file",
    help="GenBank input file. Can be compressed with gzip or bzip2"
)
parser.add_argument(
    "fasta_file", help="FASTA output datset"
)
parser.add_argument(
    "--remove_version", action="store_true",
    help="Remove version number from NCBI form formatted accession numbers. "
         "For example, this would convert 'B000657.2' to 'B000657'"
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
        print('Writing FASTA record: {}'.format(seq_id))
        print('>' + seq_id, file=output_handle)
        print(seq_record.seq, file=output_handle)
