#!/usr/bin/env python
from __future__ import print_function

import argparse
import json

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EFetch', epilog='')
    parser.add_argument('db', help='Database to use, sometimes "none" (e.g. *check)')
    parser.add_argument('dbfrom', help='Database containing input UIDs')
    parser.add_argument('cmd', choices=['neighbor', 'neighbor_score',
                                        'neighbor_history', 'acheck', 'ncheck', 'lcheck',
                                        'llinks', 'llinkslib', 'prlinks'],
                        help='ELink command mode')
    # Only used in case of neighbor_history
    #parser.add_argument('--history_out', type=argparse.FileType('w'),
    #                    help='Output history file', default='-')

    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")
    # ID Sources
    parser.add_argument('--id_xml', help='list of ids in an xml file as returned by esearch or elink')
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Fetch results from previous query')

    # Optional
    parser.add_argument('--linkname', help='Restrict results to a specific link source')
    parser.add_argument('--retmode', choices=['xml', 'json', 'uilist'], help='Output format')

    # TODO: dates, linkname, term, holding
    # neighbor or neighbor_history and dbfrom is pubmed
    # parser.add_argument('--datetype', help='Date type')
    # parser.add_argument('--reldate', help='In past N days')
    # parser.add_argument('--mindate', help='Minimum date')
    # parser.add_argument('--maxdate', help='maximum date')

    # Output
    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)
    merged_ids = c.parse_ids(args.id_list, args.id, args.history_file, args.id_xml)

    payload = {
        'dbfrom': args.dbfrom,
        'cmd': args.cmd,
    }
    if args.history_file is not None:
        payload.update(c.get_history())
    else:
        payload['id'] = ','.join(merged_ids)

    # DB can be 'none' in a few cases.
    if args.db != "none":
        payload['db'] = args.db

    if args.linkname is not None:
        payload['linkname'] = args.linkname

    #if args.retmode != "none":
    #    payload['retmode'] = args.retmode

    if (args.cmd == "neighbor_history" and args.retmode == 'json') or args.retmode == "uilist":
        payload['retmode'] = 'xml'
    elif args.retmode != "none":
        payload['retmode'] = args.retmode

    results = c.link(**payload)

    #if args.cmd == "neighbor_history":
    #    #if args.retmode == "xml":
    #    #    history = c.extract_history(results)
    #    #    args.history_out.write(json.dumps(history, indent=4))
    #    #else:
    #    #    args.history_out.write(results)

    #We're going to infer that rettype being uilist means convert to text format (which esearch does not do)
    if args.retmode is not None and args.retmode == 'uilist':
        ids = c.xmlstring2UIlist(results)
        for id in ids:
            print(id)
    elif args.cmd == "neighbor_history" and args.retmode == 'json':
        history = c.extract_history(results)
        print(json.dumps(history, indent=4))
    else:
        print(results)
