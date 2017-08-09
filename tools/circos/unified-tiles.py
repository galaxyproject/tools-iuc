#!/usr/bin/env python
import logging
import sys

from BCBio import GFF

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


# Pair up (file, extension) pairs from sys.argv
files = zip(sys.argv[1:][0::2], sys.argv[1:][1::2])


# Handlers
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

            yield (chrom, chromStart, chromEnd, lineData[4], lineData[6], lineData[5], lineData[9])


def gff3(idx, path):
    for record in GFF.parse(path):
        if len(record.features) == 0:
            continue

        for feature in sorted(record.features, key=lambda x: x.location.start):
            yield (
                record.id,
                feature.location.start,
                feature.location.end,
                feature.id or feature.qualifiers.get('Name', [None])[0],
                feature.location.strand,
                feature.qualifiers.get('score', [0.0])[0],
                feature.qualifiers.get('color', [None])[0]
            )


if __name__ == '__main__':
    for idx, (file_path, file_type) in enumerate(files):
        log.info("Processing %s.%s", file_path, file_type)

        if file_type in globals():
            func = globals()[file_type]
            for item in func(idx, file_path):
                # multiple=False
                # hs1 10292899 10301003 id=Conrad_993
                # hs1 10297766 10301003 id=Conrad_994
                lineExtra = [
                    'strand=%s' % item[4],
                    'score=%s' % item[5],
                    'value=%s' % item[5],
                ]
                if item[3] is not None:
                    lineExtra.append('id=%s' % item[3])
                if item[6] is not None:
                    lineExtra.append('color=%s' % item[6])

                sys.stdout.write(' '.join((str(item[0]), str(item[1]), str(item[2]), ','.join(lineExtra))) + '\n')
