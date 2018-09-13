#!/usr/bin/env python
import logging
import sys

import wiggle
from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


MODE = sys.argv[1]

# Pair up (file, extension) pairs from sys.argv
files = zip(sys.argv[2:][0::2], sys.argv[2:][1::2])

# Our output data structure. This could be much more efficient.
data = {}


def bed(idx, path):
    # chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    # name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    # score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
    # strand - Defines the strand - either '+' or '-'.
    # thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
    # thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
    # itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.

    with open(path, 'r') as handle:
        for line in handle:
            lineData = line.strip().split()
            chrom = lineData[0]
            chromStart = lineData[1]
            chromEnd = lineData[2]

            if chrom not in data:
                data[chrom] = {}

            for i in range(chromStart, chromEnd):
                if i not in data[chrom]:
                    data[chrom][i] = {}

                data[chrom][i][idx] = lineData[5]


# Handlers
def gff3(idx, path):
    for record in GFF.parse(path):
        if len(record.features) == 0:
            continue

        if record.id not in data:
            data[record.id] = {}

        for feature in record.features:
            if 'score' in feature.qualifiers:
                for i in range(feature.location.start, feature.location.end):
                    if i not in data[record.id]:
                        data[record.id][i] = {}

                    data[record.id][i][idx] = feature.qualifiers['score'][0]


def wig(idx, path):
    walker = wiggle.Wiggle()
    with open(path, 'r') as handle:
        for region, position, value in walker.walk(handle):
            if region not in data:
                data[region] = {}

            if position not in data[region]:
                data[region][position] = {}

            data[region][position][idx] = value


if __name__ == '__main__':
    mode_tiles_possible = True

    for idx, (file_path, file_type) in enumerate(files):
        log.info("Processing %s.%s", file_path, file_type)

        if file_type in globals():
            func = globals()[file_type]
            func(idx, file_path)

        if file_type == 'wig':
            mode_tiles_possible = False

    if MODE == 'tile' and not mode_tiles_possible:
        raise Exception("You requested a 'tile' plot with wig data, which is impossible")

    # Max number of files
    max_idx = range(len(files))

    serialized_values = None
    region_start, region_end = (None, None)

    for genome in data:
        for position in sorted(data[genome]):
            values = [
                '' if x not in data[genome][position] else data[genome][position][x]
                for x in max_idx
            ]
            if serialized_values is None:
                serialized_values = values
            if region_start is None:
                region_start = position
                region_end = position

            if values == serialized_values:
                region_end = position
            else:
                if MODE == 'histogram':
                    # histogram
                    # hs4 0 1999999 5.0000,3.0000,1.0000,19.0000
                    sys.stdout.write(' '.join(
                        (genome, str(region_start), str(region_end), ','.join(map(str, values)))
                    ) + '\n')
                elif MODE == 'heatmap':
                    # heatmap
                    # hs1 2000000 3999999 0.0000 id=hs4
                    # hs1 4000000 5999999 2.0000 id=hs1
                    # hs1 4000000 5999999 0.0000 id=hs2
                    # hs1 4000000 5999999 0.0000 id=hs3
                    # hs1 4000000 5999999 0.0000 id=hs4
                    # hs1 6000000 7999999 4.0000 id=hs2
                    for x in max_idx:
                        if x in data[genome][position]:
                            sys.stdout.write(' '.join(
                                (genome, str(region_start), str(region_end), data[genome][position][x], 'id=hm%s' % x)
                            ) + '\n')
                        else:
                            sys.stdout.write(' '.join(
                                (genome, str(region_start), str(region_end), 0.0, 'id=hm%s' % x)
                            ) + '\n')
                elif MODE == 'line':
                    # multiple=False
                    sys.stdout.write(' '.join(
                        (genome, str(region_start), str(region_end), data[genome][position][0])
                    ) + '\n')
                elif MODE == 'scatter':
                    # multiple=False
                    sys.stdout.write(' '.join(
                        (genome, str(region_start), str(region_end), data[genome][position][0])
                    ) + '\n')

                # Update start of next array
                region_start = position
                region_end = position
                # And update with new array
                serialized_values = values
