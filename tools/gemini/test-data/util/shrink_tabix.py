from __future__ import print_function

import argparse

import pysam


def main(infile, ofile, region):
    print(infile, '->', ofile)
    with pysam.Tabixfile(infile) as i:
        fformat = i.format.lower()
        if fformat == 'sam':
            fformat = 'bed'
        if ofile[-3:] == '.gz':
            ofile = ofile[:-3]
        with open(ofile, 'w') as o:
            try:
                region_it = i.fetch(region=region)
            except ValueError:
                region_it = i.fetch(region='chr' + region)
            for line in i.header:
                o.write(line + '\n')
            for line in region_it:
                o.write(str(line) + '\n')
    pysam.tabix_index(ofile, preset=fformat, force=True)


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
    args = vars(p.parse_args())
    main(**args)
