import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser()
parser.add_argument("genbank_file", help="GenBank input file")
parser.add_argument("fasta_file", help="FASTA output datset")
parser.add_argument("--remove_version", dest="remove_version", action="store_true", help="Remove version number from NCBI form formatted accessiion numbers. For example, this would convert 'B000657.2' to 'B000657'")
args = parser.parse_args()


gbk_filename = args.genbank_file
faa_filename = args.fasta_file

try:
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")
except OSError as err:
    sys.exit("OS error: {0}".format(err))

for seq_record in SeqIO.parse(input_handle, "genbank") :
    if args.remove_version:
        seq_id = seq_record.id.split('.')[0]
    else:
        seq_id = seq_record.id
    print('Writing record {}'.format( seq_id ))
    try:
        output_handle.write(">{}\n{}\n".format(seq_id, seq_record.seq))
    except OSError as err:
        sys.exit("OS error: {0}".format(err))

output_handle.close()
input_handle.close()
