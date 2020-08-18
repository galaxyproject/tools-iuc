#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import os
import sys

import eutils


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ESummary', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    # ID source
    parser.add_argument('--id_xml', help='list of ids in an xml file as returned by esearch or elink')
    parser.add_argument('--id_json', help='list of ids in a json file as returned by esearch or elink')
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Fetch results from previous query')
    parser.add_argument('--history_xml', help='Fetch results from previous query')

    # Output
    parser.add_argument('--retmode', help='Retmode')
    parser.add_argument('--retstart', type=int, default=0, help='Retstart - Starting rec number (0)')
    parser.add_argument('--retmax', type=int, default=20, help='Retmax - max number of recs returned (20, max 100000')

    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)

    payload = {
        'db': args.db,
    }

    for attr in ('retmode', 'retmax', 'retstart'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    results = []
    qkeys = []
    if args.history_file is not None or args.history_xml is not None:
        payload['retmode'] = args.retmode
        if args.history_file is not None:
            input_histories = c.get_histories()
        else:
            input_histories = c.extract_histories_from_xml_file(args.history_xml)

        for hist in input_histories:
            qkeys += [hist['query_key']]
            tmp_payload = payload
            tmp_payload.update(hist)
            results += [c.summary(**tmp_payload)]
    else:
        # There is no uilist retmode
        if args.retmode == "uilist":
            payload['retmode'] = 'xml'
        else:
            payload['retmode'] = args.retmode
        merged_ids = c.parse_ids(args.id_list, args.id, args.history_file, args.id_xml, args.id_json)
        payload['id'] = ','.join(merged_ids)
        qkeys += [1]
        results += [c.summary(**payload)]

    # There could be multiple sets of results if a history was supplied
    if args.history_file is not None or args.history_xml is not None:
        # Multiple result sets can be returned
        # Create a directory for the output files
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, r'downloads')
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)

        eprint("Writing files:")
        count = 0
        if args.retmode == 'json':
            for result in results:
                qkey = qkeys[count]
                count += 1
                file_path = os.path.join('downloads', '%s-querykey%s.json' % (args.db, qkey))
                eprint('%s-link%s.json' % (args.db, count))
                with open(file_path, 'w') as handle:
                    json_data = c.jsonstring2jsondata(result)
                    handle.write(json.dumps(json_data, indent=4))
        else:
            for result in results:
                qkey = qkeys[count]
                count += 1
                file_path = os.path.join('downloads', '%s-querykey%s.xml' % (args.db, qkey))
                eprint('%s-link%s.xml' % (args.db, count))
                with open(file_path, 'w') as handle:
                    handle.write(result)
    else:
        # When rettype is uilist, convert to text format (which elink does not do)
        if args.retmode == 'json':
            json_data = c.jsonstring2jsondata(results[0])
            print(json.dumps(json_data, indent=4))
        else:
            print(results[0])
