#!/usr/bin/env python
from __future__ import print_function

import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EInfo', epilog='')
    parser.add_argument('--db', help='Database to use')
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")
    args = parser.parse_args()

    c = eutils.Client(user_email=args.user_email, admin_email=args.admin_email)
    payload = {}
    if args.db is not None:
        payload['db'] = args.db
        payload['version'] = '2.0'
    print(c.info(**payload))
