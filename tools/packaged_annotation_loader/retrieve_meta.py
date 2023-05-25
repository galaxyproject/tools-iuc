#!/usr/bin/env python

import argparse
import json
import os

import yaml

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('galaxy_json')
    parser.add_argument(
        '-o', '--ofile',
        required=True
    )
    parser.add_argument(
        '--format', choices=['yaml', 'tab'], default='yaml'
    )
    args = parser.parse_args()

    galaxy_collection_info = json.load(open(args.galaxy_json))
    annotation_info = next(iter(galaxy_collection_info.values()))['elements']
    selected_ids = {i['name'] for i in annotation_info}
    package_meta_file = os.path.join(
        os.path.dirname(annotation_info[0]['filename']),
        'meta.yml'
    )
    meta = yaml.safe_load(open(package_meta_file))
    meta['records'] = [
        rec for rec in meta['records'] if rec['id'] in selected_ids
    ]

    with open(args.ofile, 'w') as fo:
        if args.format == 'yaml':
            yaml.dump(
                meta, fo, allow_unicode=False, default_flow_style=False
            )
        else:
            print('Annotation\tVersion', file=fo)
            for record in meta['records']:
                print(record['name'], record['version'], sep='\t', file=fo)
