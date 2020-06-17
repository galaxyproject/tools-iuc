#!/usr/bin/env python
import re
import subprocess
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature
)


def main(expterm, fasta, gff3):
    with open(fasta, 'r') as handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    # Build coords file
    with open(gff3, 'r') as handle:
        for rec in GFF.parse(handle, base_dict=seq_dict):
            with open('tmp.coords', 'w') as coords:
                for feat in rec.features:
                    if feat.type == 'gene':
                        coords.write('\t'.join([
                            feat.id,
                            str(feat.location.start + 1),
                            str(feat.location.end),
                            rec.id,
                        ]) + '\n')
            with open('tmp.fasta', 'w') as fasta_handle:
                SeqIO.write(rec, fasta_handle, 'fasta')

            cmd = ['transterm', '-p', expterm, fasta, 'tmp.coords']
            output = subprocess.check_output(cmd)
            #   TERM 1         4342 - 4366     + F    93 -11.5 -3.22878 | opp_overlap 4342, overlap 4340 4357
            ttre = re.compile(
                r'^  (?P<name>.*) (?P<start>\d+) - (?P<end>\d+)\s+'
                r'(?P<strand>[-+])\s+(?P<loc>[GFRTHNgfr]+)\s+'
                r'(?P<conf>\d+)\s+(?P<hp>[0-9.-]+)\s+(?P<tail>[0-9.-]+)'
            )

            rec.features = []
            batches = output.split('SEQUENCE ')
            for batch in batches[1:]:
                batch_lines = batch.split('\n')
                # Strip the header
                interesting = batch_lines[2:]
                unformatted = [x for x in interesting if x.startswith('  ')][0::2]
                for terminator in unformatted:
                    m = ttre.match(terminator)
                    if m:
                        start = int(m.group('start')) - 1
                        end = int(m.group('end'))
                        if m.group('strand') == '+':
                            strand = 1
                        else:
                            strand = 0

                        feature = SeqFeature(
                            FeatureLocation(start, end),
                            type="terminator",
                            strand=strand,
                            qualifiers={
                                "source": "TransTermHP_2.09",
                                "score": m.group('conf'),
                                "ID": m.group('name'),
                            }
                        )
                        rec.features.append(feature)
            yield rec


if __name__ == '__main__':
    for record in main(*sys.argv[1:4]):
        GFF.write([record], sys.stdout)
