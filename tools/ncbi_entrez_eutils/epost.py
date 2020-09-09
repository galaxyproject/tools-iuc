#!/usr/bin/env python

import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EPost', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    # ID source
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--id_json', help='list of ids in a json file as returned by esearch or elink')
    parser.add_argument('--id_xml', help='list of ids in an xml file as returned by esearch or elink')

    # Target history
    parser.add_argument('--history_xml', help='Post to new QueryKey in an existing WebEnv (XML)')
    parser.add_argument('--history_file', help='Post to new QueryKey in an existing WebEnv (JSON)')
    parser.add_argument('--webenv', help='Post to new WebEnv (History ID)')

    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)

    payload = {}
    if args.history_file is not None:
        hist = c.get_history()
        payload['WebEnv'] = hist['WebEnv']
    elif args.history_xml is not None:
        hist = c.extract_history_from_xml_file(args.history_xml)
        payload['WebEnv'] = hist['WebEnv']
    elif args.webenv is not None:
        payload['WebEnv'] = args.webenv

    merged_ids = c.parse_ids(args.id_list, args.id, None, args.id_xml, args.id_json)
    payload['id'] = ','.join(merged_ids)

    print(c.post(args.db, **payload))
