#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import sys


import eutils


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ESearch', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('term', help='Query')
    parser.add_argument('--history_file', help='Filter existing history')
    parser.add_argument('--datetype', help='Date type')
    parser.add_argument('--reldate', help='In past N days')
    parser.add_argument('--mindate', help='Minimum date')
    parser.add_argument('--maxdate', help='maximum date')
    # History
    parser.add_argument('--history_out', action="store_true", help='Output history file')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    # Output
    parser.add_argument('--retmode', help='Retmode')
    parser.add_argument('--rettype', help='Rettype')
    parser.add_argument('--retstart', type=int, default=0, help='Retstart - Starting rec number (0)')
    parser.add_argument('--retmax', type=int, default=20, help='Retmax - max number of recs returned (20, max 100000)')

    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)

    if args.retmax is not None:
        if args.retmax > 100000:
            max = 100000
        else:
            max = args.retmax

    payload = {
        'db': args.db,
        'term': args.term,
    }
    if args.history_file is not None:
        payload.update(c.get_history())

    # if args.history_out is not None:
    if args.history_out:
        payload['usehistory'] = 'y'

    payload['retmode'] = args.retmode

    for attr in ('datetype', 'reldate', 'mindate', 'maxdate', 'rettype', 'retmax', 'retstart'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    eprint("Payload used for query:", json.dumps(payload, indent=4))

    results = c.search(**payload)

    # We're going to infer that rettype being uilist means convert to text format (which esearch does not do)
    if args.retmode == 'text':
        ids = c.xmlstring2UIlist(results)
        for id in ids:
            print(id)
    elif args.retmode == 'json':
        json_data = c.jsonstring2jsondata(results)
        print(json.dumps(json_data, indent=4))
    else:
        print(results)
