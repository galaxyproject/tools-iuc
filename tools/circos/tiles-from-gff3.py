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

        for record in GFF.parse(file_path):
            if len(record.features) == 0:
                continue

            for feature in sorted(record.features, key=lambda x: x.location.start):
                # chrom chromStart chromEnd
                # name score strand
                # thickStart thickEnd itemRgb

                kv = {
                    'strand': 0 if not feature.location.strand else feature.location.strand,
                    'name': feature.qualifiers.get('Name', [None])[0] or feature.id,
                    'value': feature.qualifiers.get('score', [0])[0],
                }
                if 'color' in feature.qualifiers:
                    color = feature.qualifiers['color'][0]
                    if color not in colormap:
                        colormap[color] = 'gx-band-%s-%s' % (uniq_id, len(colormap.keys()))

                    kv['color'] = colormap[color]

                line = [
                    record.id,
                    str(int(feature.location.start)),
                    str(int(feature.location.end)),
                    ','.join(['%s=%s' % x for x in kv.items()])
                ]

                sys.stdout.write('\t'.join(line))
                sys.stdout.write('\n')
