#!/usr/bin/env python
import argparse
import eutils


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EInfo', epilog='')
    parser.add_argument('--db', help='Database to use')
    args = parser.parse_args()

    c = eutils.Client()
    payload = {}
    if args.db is not None:
        payload['db'] = args.db
    print c.info(**payload)
