#!/usr/bin/env python
import logging
import sys

from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


if __name__ == "__main__":
    attr = sys.argv[2]

    for record in GFF.parse(sys.argv[1]):
        if len(record.features) == 0:
            continue

        for feature in sorted(record.features, key=lambda x: x.location.start):
            # chrom chromStart chromEnd
            # name score strand
            # thickStart thickEnd itemRgb

            kv = {
                "strand": 0 if not feature.location.strand else feature.location.strand,
                "name": feature.qualifiers.get(attr, ["None"])[0] or feature.id,
                "value": feature.qualifiers.get("score", [0])[0],
            }

            line = [
                record.id,
                str(int(feature.location.start)),
                str(int(feature.location.end)),
                ",".join(["%s=%s" % x for x in kv.items()]),
            ]

            sys.stdout.write("\t".join(line))
            sys.stdout.write("\n")
