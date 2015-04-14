#!/usr/bin/env python
import json
import argparse
import eutils


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
    parser.add_argument('--history_out', type=argparse.FileType('w'),
                        help='Output history file')
    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file)

    payload = {
        'db': args.db,
        'term': args.term,
        'retstart': 0,
        'retmax': 20,
        # hmmm @ retmax
    }
    if args.history_file is not None:
        payload.update(c.get_history())
    if args.history_out is not None:
        payload['usehistory'] = 'y'

    for attr in ('datetype', 'reldate', 'mindate', 'maxdate'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    results = c.search(**payload)

    if args.history_out is not None:
        history = c.extract_history(results)
        args.history_out.write(json.dumps(history, indent=4))

    print results
