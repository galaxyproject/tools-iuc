#!/usr/bin/env python
import re
import os
import argparse
import subprocess
import tempfile
from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def run_reprof(query_path, modeldir):
    outtmp = tempfile.NamedTemporaryFile(delete=False)
    cmd = [
        './reprof/scripts/reprof',
        '-i', query_path,
        '--modeldir=%s' % modeldir,
        '-o', outtmp.name
    ]
    subprocess.check_call(cmd)
    outtmp.seek(0)
    data = outtmp.read()
    outtmp.close()
    os.unlink(outtmp.name)
    return data

def process_reprof_report(data):
    KEYS = ['idx', 'AA', 'PHEL', 'RI_S', 'pH', 'pE', 'pL', 'PACC', 'PREL', 'P10', 'RI_A', 'Pbe', 'Pbie']
    data_tracks = {k: [] for k in KEYS}

    for line in data.split('\n'):
        if line.startswith('#') or line.startswith('No') or len(line.strip()) == 0:
            continue

        for idx, (key, value) in enumerate(zip(KEYS, line.strip().split('\t'))):
            # numerical columns
            if idx not in (1, 2, 11, 12):
                value = int(value)

            data_tracks[key].append(value)
    return data_tracks

def storeWigData(idx, data, id, path):
    with open(path, 'a') as handle:
        handle.write('variableStep chrom=%s\n' % id)
        for (pos, val) in zip(idx, data):
            handle.write('%s %s\n' % (pos, val))

def storeGff3Data(path, id, positions, values, decodeMap):
    merged = ''.join(values)
    # http://stackoverflow.com/a/19683549
    regions = [(x[0][0].upper(), len(x[0])) for x in re.findall('((.)(\\2*))', merged)]

    location = 1

    rec = SeqRecord(Seq("ACTG"), id=id)
    for (region_char, region_length) in regions:
        # If we don't wish to decode this region, skip it.
        if region_char not in decodeMap:
            location += region_length
            continue

        region_info = decodeMap[region_char]
        # Create a feature representing this region
        region_feat = SeqFeature(
            FeatureLocation(location - 1, location - 1 + region_length),
            type=region_info['type'], strand=0,
            qualifiers={k: v for (k, v) in region_info.iteritems() if k != 'type'}
        )
        # Update our start location
        location += region_length
        rec.features.append(region_feat)

    with open(path, 'a') as handle:
        GFF.write([rec], handle)

def main(fasta, modeldir):
    for record in SeqIO.parse(fasta, 'fasta'):
        tmp = tempfile.NamedTemporaryFile(delete=False)
        SeqIO.write([record], tmp, 'fasta')
        tmp.close()

        # Run reprof
        data = process_reprof_report(run_reprof(tmp.name, modeldir))
        for col in ('RI_S', 'P10', 'RI_A', 'PACC', 'PREL', 'pH', 'pE', 'pL'):
            storeWigData(data['idx'], data[col], record.id, col + '.wig')

        storeGff3Data(
            'secondary_structure.gff3', record.id, data['idx'], data['PHEL'],
            {
                'H': {
                    'type': 'peptide_helix',
                    'label': ['Helix'],
                    'evidence': ['ECO:0000255']
                },
                'E': {
                    'type': 'beta_strand',
                    'label': ['Extended/Sheet'],
                    'evidence': ['ECO:0000255']
                },
                'L': {
                    'type': 'loop',
                    'label': ['Loop'],
                    'evidence': ['ECO:0000255']
                }
            }
        )

        storeGff3Data(
            'solvent_accessibility.gff3', record.id, data['idx'], data['Pbe'],
            {
                'B': {
                    'type': 'experimental_result_region',
                    'label': ['Buried'],
                    'evidence': ['ECO:0000255']
                },
                'E': {
                    'type': 'experimental_result_region',
                    'label': ['Exposed'],
                    'evidence': ['ECO:0000255']
                },
            }
        )


if __name__ == '__main__':
    # Grab all of the filters from our plugin loader
    parser = argparse.ArgumentParser(description='Wrapper for reprof')
    parser.add_argument('fasta', type=file, help='Fasta Input')
    parser.add_argument('modeldir')
    args = parser.parse_args()

    main(**vars(args))
