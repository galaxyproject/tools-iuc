#!/usr/bin/env python
from __future__ import print_function

import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EPost', epilog='')
    parser.add_argument('db', help='Database to use')
    parser.add_argument('--id_list', help='list of ids')
    parser.add_argument('--id', help='Comma separated individual IDs')
    parser.add_argument('--history_file', help='Post to new QueryKey in an existing WebEnv')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    args = parser.parse_args()

    c = eutils.Client(history_file=args.history_file, user_email=args.user_email, admin_email=args.admin_email)
    merged_ids = c.parse_ids(args.id_list, args.id, args.history_file)

    payload = {}
    if args.history_file is not None:
        payload.update(c.get_history())
    else:
        payload['id'] = ','.join(merged_ids)
        payload['WebEnv'] = ''

    print(c.post(args.db, **payload))
