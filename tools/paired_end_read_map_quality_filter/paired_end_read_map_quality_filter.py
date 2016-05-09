import argparse
import os
import pysam
import shutil
import tempfile


def get_tmp_file_name(prefix=None):
    tmp_dir = tempfile.mkdtemp(prefix=prefix)
    fd, tmp_name = tempfile.mkstemp(dir=tmp_dir)
    os.close(fd)
    return tmp_name

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help='Input dataset')
parser.add_argument('--output', dest='output', help='Output dataset')
args = parser.parse_args()

samfile = pysam.AlignmentFile(args.input, "rb")
tmp_output = get_tmp_file_name(prefix='mapq_filter_unsorted-')
# Samtools accepts a prefix, not a filename. It always
# adds the .bam suffix to the prefix, so we'll create a
# temporary prefix that we can use safely.
tmp_sorted_prefix = get_tmp_file_name(prefix='mapq_filter_sorted-')
pairedreads = pysam.AlignmentFile(tmp_output, "wb", template=samfile)
read1_dict = {}
read2_dict = {}

for read in samfile.fetch():
    query_name = read.query_name
    if read.is_read1:
        read1_dict.update({query_name: read})
    else:
        # read.is_read2.
        read2_dict.update({query_name: read})
    # Keep track of the Reads.
    if query_name in read1_dict and query_name in read2_dict:
        read1 = read1_dict[query_name]
        read2 = read2_dict[query_name]
        del read1_dict[query_name]
        del read2_dict[query_name]
        # Include all appropriate Reads, rescuing those that
        # would be eliminated by samtools.
        if read1.mapping_quality > 5 and read2.mapping_quality > 5:
            include = True
        elif read1.mapping_quality > 5 and read2.mapping_quality == 0:
            include = True
        elif read1.mapping_quality == 0 and read2.mapping_quality > 5:
            include = True
        else:
            include = False
        if include and read1.is_paired and not read1.is_unmapped and not read1.mate_is_unmapped and not read1.is_duplicate:
            # Save the Reads.
            pairedreads.write(read1)
            pairedreads.write(read2)

pairedreads.close()
samfile.close()
# Sort the output BAM file.
pysam.sort(tmp_output, tmp_sorted_prefix)
shutil.move('%s.bam' % tmp_sorted_prefix, args.output)
