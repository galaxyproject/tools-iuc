#!/usr/bin/env python

import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EGQuery', epilog='')
    parser.add_argument('term', help='Query')

    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")

    parser.add_argument('--version', action='version', version=eutils.Client.getVersion(), help='Version (reports Biopython version)')

    args = parser.parse_args()

    c = eutils.Client(user_email=args.user_email, admin_email=args.admin_email)

    payload = {
        'term': args.term,
    }
    results = c.gquery(**payload)
    print(results)
