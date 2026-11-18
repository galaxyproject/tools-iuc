"""Check a received MIC-50-90 reporting certificate with the program's public recipient interface.

Writes one row per declared question: the smallest and largest count the received information
allows and whether the declared answer is established. No complete count table is needed.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

from mic_50_90 import verify_reporting
from mic_50_90._version import COMPATIBLE_REPORT_VERSIONS

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('certificate', type=Path)
parser.add_argument('--table', type=Path, required=True)
parser.add_argument('--summary', type=Path, required=True)
args = parser.parse_args()

try:
    document = json.loads(args.certificate.read_text(encoding='utf-8'))
except (OSError, ValueError) as error:
    sys.exit(f'The input is not a readable JSON file: {error}')
if not isinstance(document, dict) or document.get('version') not in COMPATIBLE_REPORT_VERSIONS:
    sys.exit('The input is not a MIC-50-90 reporting certificate (reporting_certificates.json).')
records = document.get('certificates') or []
if not records:
    sys.exit('The certificate lists no collections.')

rows, lines = [], []
for record in records:
    cohort = record.get('cohort_id', 'Collection')
    result = verify_reporting(record['specification'], record['criteria'], record['disclosures'])
    lines.append(f"{cohort}: {'all declared answers are established' if result['sufficient'] else 'more information is required'}")
    for item in result['bounds']:
        rows.append({
            'cohort_id': cohort,
            'question': f"share with MIC > {item['threshold']:g} {item['unit']} {item['decision_operator']} {item['decision_fraction_text']}",
            'count_min': item['count_min'],
            'count_max': item['count_max'],
            'n': item['n'],
            'answer': item['status'],
            'collection_sufficient': 'yes' if result['sufficient'] else 'no',
        })
        lines.append(f"  MIC > {item['threshold']:g} {item['unit']}: {item['count_min']} to {item['count_max']} "
                     f"of {item['n']} isolates; {item['status']}")

with args.table.open('w', newline='', encoding='utf-8') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
lines.append('')
lines.append('The check assumes that the disclosed counts are true; it does not authenticate the source.')
args.summary.write_text('\n'.join(lines) + '\n', encoding='utf-8')
