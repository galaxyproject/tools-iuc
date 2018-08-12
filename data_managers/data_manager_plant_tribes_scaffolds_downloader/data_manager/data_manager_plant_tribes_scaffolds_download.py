#!/usr/bin/env python
#
# Data manager for downloading Plant Tribes scaffolds data.
from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tarfile
import urllib2
import zipfile


DEFAULT_DATA_TABLE_NAMES = ["plant_tribes_scaffolds"]


def add_data_table_entry(data_manager_dict, data_table_name, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get(data_table_name, [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def make_directory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def remove_directory(dir):
    if os.path.exists(dir):
        shutil.rmtree(dir)


def extract_archive(file_path, work_directory):
    if tarfile.is_tarfile(file_path):
        fh = tarfile.open(file_path, 'r:*')
    elif zipfile.is_zipfile(file_path):
        fh = zipfile.ZipFile(file_path, 'r')
    else:
        return
    fh.extractall(work_directory)


def move_files(source_directory, target_directory):
    # Move the files into defined output directory.
    for filename in os.listdir(source_directory):
        shutil.move(os.path.join(source_directory, filename), target_directory)


def url_download(url, work_directory):
    file_path = os.path.join(work_directory, os.path.basename(url))
    src = None
    dst = None
    try:
        req = urllib2.Request(url)
        src = urllib2.urlopen(req)
        dst = open(file_path, 'wb')
        while True:
            chunk = src.read(2**10)
            if chunk:
                dst.write(chunk)
            else:
                break
    except Exception as e:
        print(str(e), file=sys.stderr)
    finally:
        if src:
            src.close()
        if dst:
            dst.close()
    return file_path


def download(target_directory, web_url, config_web_url, description, data_table_names=DEFAULT_DATA_TABLE_NAMES):
    data_manager_dict = {}
    data_table_entry = {}
    # Download the scaffolds data.
    work_directory = os.path.abspath(os.path.join(os.getcwd(), 'scaffolds'))
    make_directory(work_directory)
    file_path = url_download(web_url, work_directory)
    extract_archive(file_path, work_directory)
    os.remove(file_path)
    # Move the scaffolds data files into the defined output directory.
    move_files(work_directory, target_directory)
    remove_directory(work_directory)
    # Populate the data_manager_dict with the scaffolds data entry.
    for file_path in os.listdir(target_directory):
        full_path = os.path.abspath(os.path.join(target_directory, file_path))
        entry_name = "%s" % os.path.basename(file_path)
        data_table_entry['value'] = entry_name
        data_table_entry['name'] = entry_name
        data_table_entry['path'] = full_path
        data_table_entry['description'] = description
        # Populate the data_manager_dict.
        for data_table_name in data_table_names:
            data_manager_dict = add_data_table_entry(data_manager_dict, data_table_name, data_table_entry)
    # Download the default configuration files.
    work_directory = os.path.abspath(os.path.join(os.getcwd(), 'configs'))
    make_directory(work_directory)
    file_path = url_download(config_web_url, work_directory)
    extract_archive(file_path, work_directory)
    os.remove(file_path)
    # Move the default configuration files into the defined output directory.
    source_configs_directory = os.path.join(work_directory, entry_name)
    target_configs_directory = os.path.join(target_directory, entry_name)
    move_files(source_configs_directory, target_configs_directory)
    remove_directory(work_directory)
    return data_manager_dict


parser = argparse.ArgumentParser()
parser.add_argument('--description', dest='description', default=None, help='Description')
parser.add_argument('--name', dest='name', help='Data table entry unique ID')
parser.add_argument('--out_file', dest='out_file', help='JSON output file')
parser.add_argument('--web_url', dest='web_url', help='URL for downloading scaffolds')
parser.add_argument('--config_web_url', dest='config_web_url', help='URL for downloading default configs')

args = parser.parse_args()

# Some magic happens with tools of type "manage_data" in that the output
# file contains some JSON data that allows us to define the target directory.
params = json.loads(open(args.out_file).read())
target_directory = params['output_data'][0]['extra_files_path']
make_directory(target_directory)

if args.description is None:
    description = ''
else:
    description = args.description.strip()

# Get the scaffolds data.
data_manager_dict = download(target_directory, args.web_url, args.config_web_url, description)
# Write the JSON output dataset.
with open(args.out_file, 'w') as fh:
    fh.write(json.dumps(data_manager_dict))
