# Dan Blankenberg
from __future__ import print_function

import sys


def main():
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    identifier_col = int(sys.argv[3]) - 1
    sequence_col = int(sys.argv[4]) - 1
    quality_col = int(sys.argv[5]) - 1

    max_col = max(identifier_col, sequence_col, quality_col)
    num_reads = None
    skipped_lines = 0
    with open(output_filename, 'w') as out:
        for num_reads, line in enumerate(open(input_filename)):
            fields = line.rstrip('\n\r').split('\t')
            if len(fields) > max_col:
                out.write("@%s\n%s\n+\n%s\n" % (fields[identifier_col], fields[sequence_col], fields[quality_col]))
            else:
                skipped_lines += 1

    if num_reads is None:
        print("Input was empty.")
    else:
        print("%i tabular lines were written as FASTQ reads. Be sure to use the FASTQ Groomer tool on this output before further analysis." % (num_reads + 1 - skipped_lines))


if __name__ == "__main__":
    main()
