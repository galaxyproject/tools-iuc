#!/usr/bin/env python
from __future__ import print_function

import argparse

import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ECitMatch', epilog='')
    parser.add_argument('--file', type=argparse.FileType('r'), help='Tabular file containing citations to search')

    parser.add_argument('--key', nargs='*', help='Citation Key')
    parser.add_argument('--journal_title', nargs='*', help='Journal Title')
    parser.add_argument('--year', nargs='*', help='Year')
    parser.add_argument('--volume', nargs='*', help='Volume')
    parser.add_argument('--first_page', nargs='*', help='First Page')
    parser.add_argument('--author_name', nargs='*', help='Author name')

    # Emails
    parser.add_argument('--user_email', help="User email")
    parser.add_argument('--admin_email', help="Admin email")
    args = parser.parse_args()

    c = eutils.Client(user_email=args.user_email, admin_email=args.admin_email)

    citations = []
    if args.file is None:
        for key, journal, year, volume, first_page, author_name in \
                zip(args.key, args.journal_title, args.year, args.volume, args.first_page, args.author_name):
            citations.append({
                'key': key,
                'journal': journal,
                'year': year,
                'volume': volume,
                'first_page': first_page,
                'author_name': author_name,
            })
    else:
        for line in args.file:
            line = line.strip()
            if not line.startswith('#'):
                tmp = line.split('\t')
                try:
                    citations.append({
                        'journal': tmp[0],
                        'year': tmp[1],
                        'volume': tmp[2],
                        'first_page': tmp[3],
                        'author_name': tmp[4],
                        'key': tmp[5],
                    })
                except KeyError:
                    print("Could not parse line: %s" % line)

    payload = {
        'db': 'pubmed',
        'bdata': citations
    }

    results = c.citmatch(**payload)
    # We get data back as pipe separated, so just replace those with tabs
    print(results.replace('|', '\t'))
