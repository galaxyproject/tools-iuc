#!/usr/bin/env python
from __future__ import print_function

import argparse
import json
import os.path
import subprocess
import tarfile
import tempfile
import zipfile
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import urlopen


def url_download(url, workdir):
    file_path = os.path.join(workdir, 'download.dat')
    src = None
    dst = None
    try:
        src = urlopen(url)
        with open(file_path, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    finally:
        if src:
            src.close()
    if tarfile.is_tarfile(file_path):
        fh = tarfile.open(file_path, 'r:*')
    elif zipfile.is_zipfile(file_path):
        fh = zipfile.ZipFile(file_path, 'r')
    else:
        return
    fh.extractall(workdir)
    os.remove(file_path)


def cat_prepare(install_dir):
    cmd = ['CAT', 'prepare' '--fresh', '-q']
    cmd_stdout = tempfile.NamedTemporaryFile()
    cmd_stderr = tempfile.NamedTemporaryFile()
    return_code = subprocess.call(cmd, shell=True, cwd=install_dir,
                                  stdout=cmd_stdout, stderr=cmd_stderr)
    if return_code:
        msg = "stdout:\n%s\nstderr:\n%s" % (cmd_stdout.read(),
                                            cmd_stderr.read())
        cmd_stdout.close()
        cmd_stderr.close()
        raise Exception('Error: (%s), returncode=%s %s'
                        % (' '.join(cmd), return_code, msg))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file')
    parser.add_argument('--install_path')
    parser.add_argument('--db_url', default=None)
    args = parser.parse_args()

    if not os.path.exists(args.install_path):
        os.makedirs(args.install_path)
    if args.db_url:
        url_download(args.db_url, args.install_path)
    else:
        cat_prepare(args.install_path)

    cat_path = None
    cat_db = None
    tax_db = None
    for root, dirs, files in os.walk(args.install_path):
        for dname in dirs:
            if dname.endswith('CAT_database'):
                cat_db = dname
            elif dname.endswith('taxonomy'):
                tax_db = dname
        if cat_db and tax_db:
            cat_path = root
            break

    cat_dir = os.path.basename(cat_path)
    # params = json.loads(open(args.config_file).read())
    dm_dict = {}
    dm_dict['data_tables'] = dm_dict.get('data_tables', {})
    data_table = 'cat_database'
    dm_dict['data_tables'][data_table]\
        = dm_dict['data_tables'].get(data_table, [])
    data_table_entry = dict(value=cat_dir, name=cat_dir,
                            database_folder=os.path.join(cat_dir, cat_db),
                            taxonomy_folder=os.path.join(cat_dir, tax_db))
    dm_dict['data_tables'][data_table].append(data_table_entry)
    # save info to json file
    open(args.config_file, 'wb').write(json.dumps(dm_dict))


if __name__ == "__main__":
    main()
