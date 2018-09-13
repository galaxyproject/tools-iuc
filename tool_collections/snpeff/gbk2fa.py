import argparse
import bz2
import contextlib
import gzip
import sys

import magic
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("genbank_file", help="GenBank input file. Can be compressed with gzip or bzip2")
parser.add_argument("fasta_file", help="FASTA output datset")
parser.add_argument("--remove_version", dest="remove_version", action="store_true", help="Remove version number from NCBI form formatted accession numbers. For example, this would convert 'B000657.2' to 'B000657'")
args = parser.parse_args()

gbk_filename = args.genbank_file
fa_filename = args.fasta_file


@contextlib.contextmanager
def get_file_handle(gbk_filename):
    f_type = magic.from_file(args.genbank_file, mime=True)
    if f_type == 'text/plain':
        input_handle = open(gbk_filename, "r")
    elif f_type == 'application/gzip' or f_type == 'application/x-gzip':
        input_handle = gzip.open(gbk_filename, "rt")
    elif f_type == 'application/x-bzip2':
        input_handle = bz2.open(gbk_filename, "rt")
    else:
        sys.exit("Cannot process file of type {}. Only plain, gzip'ed, and bzip2'ed genbank files are accepted ".format(f_type))
    yield input_handle
    input_handle.close()


with get_file_handle(gbk_filename) as input_handle, open(fa_filename, "w") as output_handle:

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        if args.remove_version:
            seq_id = seq_record.id.split('.')[0]
        else:
            seq_id = seq_record.id
        print('Writing FASTA record: {}'.format( seq_id ))
        output_handle.write(">{}\n{}\n".format(seq_id, seq_record.seq))
