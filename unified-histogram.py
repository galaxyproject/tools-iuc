#!/usr/bin/env python
import logging
import sys
logging.basicConfig(level=logging.INFO)
log = logging.getLogger()
from BCBio import GFF

# Pair up (file, extension) pairs from sys.argv
files = zip(sys.argv[1:][0::2], sys.argv[1:][1::2])

# Our output data structure. This could be much more efficient.
data = {}

# Handlers
def gff3(idx, path):
    for record in GFF.parse(path):
        if len(record.features) == 0:
            continue

        if record.id not in data:
            data[record.id] = {}

        for feature in record.features:
            if 'score' in feature.qualifiers:
                for i in xrange(feature.location.start, feature.location.end):
                    if i not in data[record.id]:
                        data[record.id][i] = {}

                    data[record.id][i][idx] = feature.qualifiers['score'][0]

for idx, (file_path, file_type) in enumerate(files):
    log.info("Processing %s.%s", file_path, file_type)
    if file_type in globals():
        func = globals()[file_type]
        func(idx, file_path)

# Max number of files
max_idx = range(len(files))

serialized_values = None
region_start, region_end = (None, None)

for genome in data:
    for position in sorted(data[genome]):
        values = ','.join([
            '' if x not in data[genome][position] else data[genome][position][x]
            for x in max_idx
        ])
        if serialized_values is None:
            serialized_values = values
        if region_start is None:
            region_start = position
            region_end = position

        if values == serialized_values:
            region_end = position
        else:
            print genome, region_start, region_end + 1, values
            # Update start of next array
            region_start = position
            region_end = position
            # And update with new array
            serialized_values = values
        # hs4 0 1999999 5.0000,3.0000,1.0000,19.0000
