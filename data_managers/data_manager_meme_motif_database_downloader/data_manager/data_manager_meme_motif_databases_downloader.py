#!/usr/bin/env python
#
# Data manager for downloading MEME Motif Databases.
import argparse
import json
import os
import shutil
import sys
import tarfile
import urllib2
import zipfile

DEFAULT_DATA_TABLE_NAMES = ["meme_motif_databases"]

parser = argparse.ArgumentParser()
parser.add_argument('--description', dest='description', default=None, help='Description')
parser.add_argument('--name', dest='name', help='Data table entry unique ID')
parser.add_argument('--out_file', dest='out_file', help='JSON output file')
parser.add_argument('--web_url', dest='web_url', help='URL for downloading MEME motif databases')

args = parser.parse_args()


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
        sys.stderr.write('An error occurred: %s' % e)
    finally:
        if src:
            src.close()
        if dst:
            dst.close()
    return file_path


def download(target_directory, web_url, name, description, data_table_names=DEFAULT_DATA_TABLE_NAMES):
    motif_databases_directory = os.path.join(target_directory, 'motif_databases')
    data_manager_dict = {}
    # Download the databases.
    work_directory = os.path.abspath(os.path.join(os.getcwd(), 'meme_motif_databases'))
    make_directory(work_directory)
    file_path = url_download(web_url, work_directory)
    extract_archive(file_path, work_directory)
    os.remove(file_path)
    # Move the database files into the defined output directory.
    move_files(work_directory, target_directory)
    remove_directory(work_directory)
    # Populate the data_manager_dict with the database data entries.
    for file_name in os.listdir(motif_databases_directory):
        data_table_entry = {}
        full_path = os.path.abspath(os.path.join(motif_databases_directory, file_name))
        # Eliminate anything that is not a directory.
        if os.path.isdir(full_path):
            entry_name = "%s" % os.path.basename(full_path)
            data_table_entry['value'] = entry_name
            data_table_entry['name'] = entry_name
            data_table_entry['path'] = full_path
            data_table_entry['description'] = description
            # Populate the data_manager_dict.
            for data_table_name in data_table_names:
                data_manager_dict = add_data_table_entry(data_manager_dict, data_table_name, data_table_entry)
    return data_manager_dict

params = json.loads(open(args.out_file).read())
target_directory = params['output_data'][0]['extra_files_path']
make_directory(target_directory)

if args.description is None:
    description = ''
else:
    description = args.description.strip()

# Get the databases.
data_manager_dict = download(target_directory, args.web_url, args.name, description)
# Write the JSON output dataset.
fh = open(args.out_file, 'wb')
fh.write(json.dumps(data_manager_dict))
fh.close()
