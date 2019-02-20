#!/usr/bin/env python
import argparse
import logging
import sys

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature
)
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def parse_xmfa(xmfa):
    """Simple XMFA parser until https://github.com/biopython/biopython/pull/544
    """
    current_lcb = []
    current_seq = {}
    for line in xmfa.readlines():
        if line.startswith('#'):
            continue

        if line.strip() == '=':
            if 'id' in current_seq:
                current_lcb.append(current_seq)
                current_seq = {}
            yield current_lcb
            current_lcb = []
        else:
            line = line.strip()
            if line.startswith('>'):
                if 'id' in current_seq:
                    current_lcb.append(current_seq)
                    current_seq = {}
                data = line.strip().split()
                id, loc = data[1].split(':')
                start, end = loc.split('-')
                current_seq = {
                    'rid': '_'.join(data[1:]),
                    'id': id,
                    'start': int(start),
                    'end': int(end),
                    'strand': 1 if data[2] == '+' else -1,
                    'seq': ''
                }
            else:
                current_seq['seq'] += line.strip()


def _percent_identity(a, b):
    """Calculate % identity, ignoring gaps in the host sequence
    """
    match = 0
    mismatch = 0
    for char_a, char_b in zip(list(a), list(b)):
        if char_a == '-':
            continue
        if char_a == char_b:
            match += 1
        else:
            mismatch += 1

    if match + mismatch == 0:
        return 0
    return 100 * float(match) / (match + mismatch)


def _id_tn_dict(sequences):
    """Figure out sequence IDs
    """
    label_convert = {}
    if sequences is not None:
        if len(sequences) == 1:
            for i, record in enumerate(SeqIO.parse(sequences[0], 'fasta')):
                label_convert[str(i + 1)] = record.id
        else:
            for i, sequence in enumerate(sequences):
                for record in SeqIO.parse(sequence, 'fasta'):
                    label_convert[str(i + 1)] = record.id
                    continue
    return label_convert


def convert_xmfa_to_gff3(xmfa_file, relative_to='1', sequences=None, window_size=1000):
    label_convert = _id_tn_dict(sequences)

    lcbs = parse_xmfa(xmfa_file)

    records = [SeqRecord(Seq("A"), id=label_convert.get(relative_to, relative_to))]
    for lcb in lcbs:
        ids = [seq['id'] for seq in lcb]

        # Doesn't match part of our sequence
        if relative_to not in ids:
            continue

        # Skip sequences that are JUST our "relative_to" genome
        if len(ids) == 1:
            continue

        parent = [seq for seq in lcb if seq['id'] == relative_to][0]
        others = [seq for seq in lcb if seq['id'] != relative_to]

        for other in others:
            targetlist = [label_convert.get(other['id'], other['id']), str(other['start']), str(other['end'])]
            if other['strand'] != '.':
                targetlist.append(str(other['strand']))
            other['feature'] = SeqFeature(
                FeatureLocation(parent['start'], parent['end'] + 1),
                type="match", strand=parent['strand'],
                qualifiers={
                    "source": "progressiveMauve",
                    "Target": " ".join(targetlist),
                    "ID": label_convert.get(other['id'], 'xmfa_' + other['rid'])
                }
            )

        for i in range(0, len(lcb[0]['seq']), window_size):
            block_seq = parent['seq'][i:i + window_size]
            real_window_size = len(block_seq)
            real_start = abs(parent['start']) - parent['seq'][0:i].count('-') + i
            real_end = real_start + real_window_size - block_seq.count('-')

            if (real_end - real_start) < 10:
                continue

            if parent['start'] < 0:
                strand = -1
            else:
                strand = 1

            for other in others:
                pid = _percent_identity(block_seq, other['seq'][i:i + real_window_size])
                # Ignore 0% identity sequences
                if pid == 0:
                    continue

                # Support for Biopython 1.68 and above, which removed sub_features
                if not hasattr(other['feature'], "sub_features"):
                    other['feature'].sub_features = []
                other['feature'].sub_features.append(
                    SeqFeature(
                        FeatureLocation(real_start, real_end),
                        type="match_part", strand=strand,
                        qualifiers={
                            "source": "progressiveMauve",
                            'score': pid
                        }
                    )
                )

        for other in others:
            records[0].features.append(other['feature'])
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert XMFA alignments to gff3', prog='xmfa2gff3')
    parser.add_argument('xmfa_file', type=file, help='XMFA File')
    parser.add_argument('--window_size', type=int, help='Window size for analysis', default=1000)
    parser.add_argument('--relative_to', type=str, help='Index of the parent sequence in the MSA', default='1')
    parser.add_argument('--sequences', type=file, nargs='+',
                        help='Fasta files (in same order) passed to parent for reconstructing proper IDs')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    result = convert_xmfa_to_gff3(**vars(args))
    GFF.write(result, sys.stdout)
