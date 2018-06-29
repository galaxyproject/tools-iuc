#!/usr/bin/env python
import logging
import sys
import hashlib

from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == '__main__':
    colormap = {}
    # We need to ensure the colour maps do not overlap (in name) across
    # multiple runs of this script. That would cause a circos error.
    # The best, reproducible idea I've had so far is hashing.
    input_str = '|'.join(sys.argv[1:])
    uniq_id = hashlib.md5(str(input_str).encode('utf-8')).hexdigest()[0:6]

    for idx, file_path in enumerate(sys.argv[1:]):
        log.info("Processing %s", file_path)

        with open(file_path, 'r') as handle:
            for line in handle:
                lineData = line.strip().split()
                # BED3+ chrom chromStart chromEnd
                # BED6+ name score strand
                # BED9+ thickStart thickEnd itemRgb
                kv = {}
                if len(lineData) >= 9:
                    color = lineData[8]
                    if color not in colormap:
                        colormap[color] = 'gx-band-%s-%s' % (uniq_id, len(colormap.keys()))

                    kv['color'] = colormap[color]
                if len(lineData) >= 6:
                    kv['strand'] = lineData[5].replace('+', '1').replace('-', '-1')
                    kv['name'] = lineData[3]
                    kv['value'] = lineData[4]

                line = [
                    lineData[0],  # chrom
                    lineData[1],  # chromStart
                    lineData[2],  # chromEnd
                    ','.join(['%s=%s' % x for x in kv.items()])
                ]

                sys.stdout.write('\t'.join(line))
                sys.stdout.write('\n')
