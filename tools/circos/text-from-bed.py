#!/usr/bin/env python
import logging
import sys

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as handle:
        for line in handle:
            lineData = line.strip().split()
            # BED3+ chrom chromStart chromEnd
            # BED6+ name score strand
            # BED9+ thickStart thickEnd itemRgb
            kv = {}
            if len(lineData) >= 6:
                kv['strand'] = lineData[5].replace('+', '1').replace('-', '-1')
                kv['value'] = lineData[4]
            else:
                sys.exit("Must be BED6+")

            line = [
                lineData[0],  # chrom
                lineData[1],  # chromStart
                lineData[2],  # chromEnd
                lineData[3],
                ','.join(['%s=%s' % x for x in kv.items()])
            ]

            sys.stdout.write('\t'.join(line))
            sys.stdout.write('\n')
