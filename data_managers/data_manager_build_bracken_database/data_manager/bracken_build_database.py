#!/usr/bin/env python

import argparse
import errno
import json
import os
import subprocess
import uuid
from sys import stderr, stdout


DATA_TABLE_NAME = "bracken_databases"


def bracken_build_database(target_directory, bracken_build_args, database_name, prebuilt=False, data_table_name=DATA_TABLE_NAME):

    database_value = str(uuid.uuid4())

    database_kmer_name = 'database' + str(bracken_build_args['read_len']) + 'mers.kmer_distrib'

    database_path_name = bracken_build_args['kraken_database']

    database_path = target_directory

    database_kmer_path = os.path.join(target_directory, database_kmer_name)

    if not prebuilt:
        bracken_build_args_list = [
            '-t', bracken_build_args['threads'],
            '-k', bracken_build_args['kmer_len'],
            '-l', bracken_build_args['read_len'],
            '-d', database_path
        ]

        cp = subprocess.run(['bracken-build'] + bracken_build_args_list,
                            capture_output=True, text=True)
        if cp.stdout:
            stdout.write(cp.stdout)
        if cp.stderr:
            stderr.write(cp.stderr)
        cp.check_returncode()
        if cp.stdout and 'ERROR:' in cp.stdout:
            raise Exception("Error in calling bracken-build, but returncode is '%i'." % cp.returncode)
    else:
        if not os.path.exists(database_kmer_path):
            raise Exception("Requested prebuilt database does not include .kmer_distrib file: '%s'" % database_kmer_path)

    data_table_entry = {
        "data_tables": {
            data_table_name: [
                {
                    "value": database_value,
                    "name": database_name,
                    "path": database_kmer_name,
                    "db_dir": database_path_name
                }
            ]
        }
    }

    return data_table_entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument('--threads', dest='threads', default=1, help='threads')
    parser.add_argument('--kmer-len', dest='kmer_len', help='K-mer length')
    parser.add_argument('--read-len', dest='read_len', help='Read length')
    parser.add_argument('--kraken-db', dest='kraken_database', help='Kraken Database')
    parser.add_argument('--database-name', dest='database_name', help='Database Name')
    parser.add_argument('--root-dir', dest='root_dir', help='Root directory for building')
    parser.add_argument('--prebuilt', action='store_true', dest='prebuilt', help='Use pre-built DB')
    args = parser.parse_args()

    if args.prebuilt:
        bracken_build_args = {
            'threads': args.threads,
            'read_len': args.read_len,
            'kraken_database': args.kraken_database,
        }
    else:
        bracken_build_args = {
            'threads': args.threads,
            'kmer_len': args.kmer_len,
            'read_len': args.read_len,
            'kraken_database': args.kraken_database,
        }

    try:
        os.mkdir(args.root_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(args.root_dir):
            pass
        else:
            raise

    data_manager_output = {}

    data_manager_output = bracken_build_database(
        args.root_dir,
        bracken_build_args,
        args.database_name,
        args.prebuilt,
    )

    with open(args.data_manager_json, 'w') as fh:
        json.dump(data_manager_output, fh, sort_keys=True)


if __name__ == "__main__":
    main()
