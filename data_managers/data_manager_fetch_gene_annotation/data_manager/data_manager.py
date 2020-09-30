#!/usr/bin/env python

import argparse
import json
import os
import sys
from urllib.request import Request, urlopen


def url_download(url, workdir):
    # Attempt to download gene annotation file from a given url
    file_path = os.path.abspath(os.path.join(workdir, os.path.basename(url)))
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(file_path, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src:
            src.close()
    return file_path


def download(value, name, description, url, out_file):
    with open(out_file) as fh:
        params = json.loads(fh.read())
    workdir = params['output_data'][0]['extra_files_path']
    os.makedirs(workdir)
    file_path = url_download(url, workdir)
    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = value
    data_manager_entry['name'] = name
    data_manager_entry['description'] = description
    data_manager_entry['path'] = file_path
    data_manager_json["data_tables"]["gff_gene_annotations"] = data_manager_entry
    with open(out_file, 'w') as fh:
        fh.write(json.dumps(data_manager_json, sort_keys=True))


parser = argparse.ArgumentParser()
parser.add_argument('--value', dest='value', action='store', help='Data table entry unique ID')
parser.add_argument('--description', dest='description', action='store', help='Description entry')
parser.add_argument('--name', dest='name', action='store', help='Name')
parser.add_argument('--out_file', dest='out_file', action='store', help='JSON filename')
parser.add_argument('--url', dest='url', action='store', help='Url to download gtf file from')

args = parser.parse_args()

download(args.value, args.name, args.description, args.url, args.out_file)
