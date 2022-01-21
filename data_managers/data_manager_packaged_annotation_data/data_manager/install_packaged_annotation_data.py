#!/usr/bin/env python

import argparse
import datetime
import json
import os
import re
from urllib.request import urlretrieve

import yaml


class PackagedAnnotationMeta():
    @classmethod
    def from_file(cls, fname):
        meta = yaml.safe_load(open(fname))
        return cls(meta)

    def __init__(self, meta_dict):
        if 'build' not in meta_dict:
            meta_dict['build'] = datetime.date.today().isoformat()
        if 'volume' not in meta_dict:
            meta_dict['volume'] = 1

        required_meta = ['name', 'build', 'volume', 'refgenome', 'records']
        for key in required_meta:
            if not meta_dict.get(key):
                raise KeyError(
                    'Required info "{0}" missing from metadata'
                    .format(key)
                )
        required_record_meta = ['id', 'name', 'version', 'format', 'source']
        for key in required_record_meta:
            for record in meta_dict['records']:
                if not record.get(key):
                    raise KeyError(
                        '{0}\n'
                        'Required info "{1}" missing from record metadata'
                        .format(record, key)
                    )
        self.meta = meta_dict
        self.meta['id'] = self._get_id()

    def _get_id(self):
        components = [
            self.meta['name'],
            self.meta['refgenome'],
            str(self.meta['volume']),
            str(self.meta['build'])
        ]
        return '__'.join(
            [
                re.sub(r'[^a-zA-Z_0-9\-]', '', i.replace(' ', '_'))
                for i in components
            ]
        )

    def records(self, full_record_names=False):
        for record in self.meta['records']:
            ret = record.copy()
            if full_record_names:
                ret['name'] = self._full_record_name(record)
            yield ret

    def fullname(self):
        return '{0} ({1}, vol:{2}/build:{3})'.format(
            self.meta['name'],
            self.meta['refgenome'],
            self.meta['volume'],
            self.meta['build']
        )

    def _full_record_name(self, record):
        return '{0} ({1}, {2}; from {3}/vol:{4}/build:{5})'.format(
            record['name'], record['version'],
            self.meta['refgenome'],
            self.meta['name'],
            self.meta['volume'],
            self.meta['build']
        )

    def dump(self, fname):
        with open(fname, 'w') as fo:
            yaml.dump(
                self.meta, fo, allow_unicode=False, default_flow_style=False
            )


def fetch_data(source_url, target_file):
    final_file, headers = urlretrieve(source_url, target_file)


def meta_to_dm_records(meta, dbkey=None):
    data_table_rows = []
    for record in meta.records(full_record_names=True):
        data_table_rows.append(
            {
                'value': '{0}:{1}'.format(meta.meta['id'], record['id']),
                'dbkey': dbkey or meta.meta['refgenome'],
                'data_name': record['name'],
                'data_id': record['id'],
                'data_format': record['format'],
                'package_id': meta.meta['id'],
                'package_name': meta.fullname(),
                'path': '{0}/{1}'.format(
                    meta.meta['volume'],
                    meta.meta['build']
                )
            }
        )
    return data_table_rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata')
    parser.add_argument(
        '-o', '--galaxy-datamanager-json',
        required=True
    )
    parser.add_argument('-t', '--target-directory', default=None)
    parser.add_argument('--dbkey', default=None)
    args = parser.parse_args()

    if args.target_directory:
        if not os.path.isdir(args.target_directory):
            os.mkdir(args.target_directory)
    else:
        args.target_directory = os.getcwd()

    meta = PackagedAnnotationMeta.from_file(args.metadata)

    for record in meta.records():
        fetch_data(
            record['source'],
            os.path.join(args.target_directory, record['id'])
        )

    meta.dump(os.path.join(args.target_directory, 'meta.yml'))

    # Finally, we prepare the metadata for the new data table record ...
    data_manager_dict = {
        'data_tables': {
            'packaged_annotation_data': meta_to_dm_records(meta, args.dbkey)
        }
    }

    # ... and save it to the json results file
    with open(args.galaxy_datamanager_json, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)
