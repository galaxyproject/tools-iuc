#!/usr/bin/env python
import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EFetch', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    # ID source
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Fetch results from previous query')

    # Output
    parser.add_argument('--retmode', help='Retmode')
    parser.add_argument('--rettype', help='Rettype')
    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)
    merged_ids = c.parse_ids(args.id_list, args.id, args.history_file)

    payload = {}
    if args.history_file is not None:
        payload.update(c.get_history())
    else:
        payload['id'] = ','.join(merged_ids)

    for attr in ('retmode', 'rettype'):
        if getattr(args, attr, None) is not None:
            payload[attr] = getattr(args, attr)

    c.fetch(args.db, ftype=args.retmode, **payload)
