#!/usr/bin/env python

import argparse
import json
import logging
import os

import eutils


logging.basicConfig(level=logging.DEBUG)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EFetch', epilog='')
    parser.add_argument('db', help='Database to use, sometimes "none" (e.g. *check)')
    parser.add_argument('dbfrom', help='Database containing input UIDs')
    parser.add_argument('cmd', choices=['neighbor', 'neighbor_score',
                                        'neighbor_history', 'acheck', 'ncheck', 'lcheck',
                                        'llinks', 'llinkslib', 'prlinks'],
                        help='ELink command mode')

    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    # ID Sources
    parser.add_argument('--id_xml', help='list of ids in an xml file as returned by esearch or elink')
    parser.add_argument('--id_json', help='list of ids in a json file as returned by esearch or elink')
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Fetch results from previous query')
    parser.add_argument('--history_xml', help='Fetch results from previous query')

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

    payload = {
        'dbfrom': args.dbfrom,
        'cmd': args.cmd,
    }

    # DB can be 'none' in a few cases.
    if args.db != "none":
        payload['db'] = args.db

    if args.linkname is not None:
        payload['linkname'] = args.linkname

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
            results += [c.link(**tmp_payload)]
    else:
        # There is no uilist retmode
        if args.retmode == "uilist":
            payload['retmode'] = 'xml'
        else:
            payload['retmode'] = args.retmode
        merged_ids = c.parse_ids(args.id_list, args.id, args.history_file, args.id_xml, args.id_json)
        payload['id'] = ','.join(merged_ids)
        qkeys += [1]
        results += [c.link(**payload)]

    # There could be multiple sets of results if a history was supplied
    if args.history_file is not None or args.history_xml is not None:
        # Multiple result sets can be returned
        # Create a directory for the output files
        current_directory = os.getcwd()
        final_directory = os.path.join(current_directory, r'downloads')
        if not os.path.exists(final_directory):
            os.makedirs(final_directory)

        logging.info("Writing files:")
        # When rettype is uilist, convert to text format (which elink does not do)
        count = 0
        if args.retmode == 'uilist':
            for result in results:
                qkey = qkeys[count]
                count += 1
                ids = c.xmlstring2UIlist(result)
                file_path = os.path.join('downloads', '%s-querykey%s.tabular' % (args.db, qkey))
                logging.info('%s.tabular' % (args.db))
                with open(file_path, 'w') as handle:
                    for id in ids:
                        handle.write(id)
                        handle.write(os.linesep)
        elif args.retmode == 'json':
            for result in results:
                qkey = qkeys[count]
                count += 1
                file_path = os.path.join('downloads', '%s-querykey%s.json' % (args.db, qkey))
                logging.info('%s-link%s.json' % (args.db, count))
                with open(file_path, 'w') as handle:
                    json_data = c.jsonstring2jsondata(result)
                    handle.write(json.dumps(json_data, indent=4))
        else:
            for result in results:
                qkey = qkeys[count]
                count += 1
                file_path = os.path.join('downloads', '%s-querykey%s.xml' % (args.db, qkey))
                logging.info('%s-link%s.xml' % (args.db, count))
                with open(file_path, 'w') as handle:
                    handle.write(result)
    else:
        # When rettype is uilist, convert to text format (which elink does not do)
        if args.retmode == 'uilist':
            ids = c.xmlstring2UIlist(results[0])
            for id in ids:
                print(id)
        elif args.retmode == 'json':
            json_data = c.jsonstring2jsondata(results[0])
            print(json.dumps(json_data, indent=4))
        else:
            print(results[0])
