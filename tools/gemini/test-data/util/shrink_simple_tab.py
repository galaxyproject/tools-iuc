from __future__ import print_function

import argparse

from functools import partial


def keep_line(line, pos_cols, region):
    fields = line.rstrip().split(b'\t')
    if fields[pos_cols[0]] == region[0]:  # same chromosome
        if (
            region[1] < int(fields[pos_cols[1]]) < region[2]
        ) or (
            region[1] < int(fields[pos_cols[2]]) < region[2]
        ):
            return True


def main(infile, ofile, num_header_lines):
    print(infile, '->', ofile)
    with open(infile, 'rb') as i:
        with open(ofile, 'wb') as o:
            # copy header lines
            for c in range(num_header_lines):
                o.write(next(i))
            for line in i:
                if keep_line(line):
                    o.write(line)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('infile')
    p.add_argument(
        '-r', '--region',
        required=True,
        help='the region of the input file to rewrite'
    )
    p.add_argument(
        '-o', '--ofile',
        required=True,
        help="the name of the output file"
    )
    p.add_argument(
        '-c', '--cols',
        nargs=3, type=int, required=True,
        help="the columns of the input file specifying chrom, start and stop, "
             "respectively"
    )
    p.add_argument(
        '-n', '--num-header-lines',
        type=int, default=0,
        help='the number of header lines present in the input; These will '
             'always be copied over to the new file.'
    )
    args = vars(p.parse_args())

    chrom, reg = args['region'].split(':')
    region = [chrom.encode()] + [int(x) for x in reg.split('-')]
    keep_line = partial(keep_line, pos_cols=args['cols'], region=region)

    main(args['infile'], args['ofile'], args['num_header_lines'])
