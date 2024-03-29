#!/usr/bin/env python
# painfully converted from b0rken perl from
# https://unpkg.com/browse/jbrowse-plugin-mafviewer@1.0.6/dist/
# license is Apache2_license.txt included here

import sys

id = 0
buffer = ''
start = 0
end = 0
score = 0
chrom = ''

db = "%s." % sys.argv[1]
# Read input from stdin
for line in sys.stdin:
    line = line.strip()
    if not line:
        continue

    line = line.split()
    if line[0] == 's' and line[1].startswith(db):
        chrom = line[1]
        chrom = chrom.replace(db, '')
        start = int(line[2])
        end = int(line[2]) + int(line[3])
        line = line[1:]
        line = ':'.join(line)
        temp = line
        buffer = temp if buffer == '' else f"{buffer},{temp}"
    elif line[0] == 'a':
        score = int(line[1].split('=')[1])
        if id > 0:
            sys.stdout.write('\t'.join([chrom, '%d' % start, '%d' % end, f"{sys.argv[1]}_{id}", '%d' % score, buffer]) + '\n')
        id += 1
        buffer = ''
    elif line[0] == 's':
        line = line[1:]
        line = ':'.join(line)
        temp = line
        buffer = temp if buffer == '' else f"{buffer},{temp}"

sys.stdout.write('\t'.join([chrom, '%d' % start, '%d' % end, f"{sys.argv[1]}_{id}", '%d' % score, buffer]) + '\n')
