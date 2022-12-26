#!/usr/bin/env python

import argparse
import json
import os
import shutil
import sys
import tarfile
from urllib.parse import urlparse
from urllib.request import Request
from urllib.request import urlopen


def url_download(url, target_directory):
    url_parts = urlparse(url)
    tarball = os.path.abspath(os.path.join(target_directory, os.path.basename(url_parts.path)))
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(tarball, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    except Exception as e:
        sys.exit(str(e))
    finally:
        if src is not None:
            src.close()
    if tarfile.is_tarfile(tarball):
        fh = tarfile.open(tarball, 'r:*')
    else:
        return tarball
    fh.extractall(target_directory)
    fh.close()
    os.remove(tarball)
    # The tarball extraction will create a directory named
    # something like release202 in the target_directory, so
    # we need to move the items in that directory to the
    # target directory.
    subdir = next(os.walk(target_directory))[1][0]
    subdir_path = os.path.join(target_directory, subdir)
    items = os.listdir(subdir_path)
    for item in items:
        item_path = os.path.join(subdir_path, item)
        shutil.move(item_path, target_directory)
    os.rmdir(subdir_path)
    return target_directory


def download(database_id, database_name, url, out_file):

    with open(out_file) as fh:
        params = json.load(fh)

    target_directory = params['output_data'][0]['extra_files_path']
    os.makedirs(target_directory)
    file_path = url_download(url, target_directory)

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = database_id
    data_manager_entry['name'] = database_name
    data_manager_entry['path'] = file_path
    data_manager_json["data_tables"]["gtdbtk_database"] = data_manager_entry

    with open(out_file, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


parser = argparse.ArgumentParser()

parser.add_argument('--database_name', dest='database_name', help='GTDB-Tk database display name')
parser.add_argument('--database_id', dest='database_id', help='Unique GTDB-Tk database id')
parser.add_argument('--url', dest='url', help='URL to download GTDB-Tk databse version')
parser.add_argument('--out_file', dest='out_file', help='JSON output file')

args = parser.parse_args()

download(args.database_id, args.database_name, args.url, args.out_file)
