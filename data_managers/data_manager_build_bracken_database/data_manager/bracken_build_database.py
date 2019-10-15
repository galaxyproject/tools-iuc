#!/usr/bin/env python

from __future__ import print_function

import argparse
import errno
import json
import os
import subprocess
import uuid


DATA_TABLE_NAME = "bracken_databases"


def bracken_build_database(target_directory, bracken_build_args, database_name, data_table_name=DATA_TABLE_NAME):

    database_value = str(uuid.uuid4())

    database_name = database_name

    database_path = os.path.join(bracken_build_args['kraken_database'], 'database' + str(bracken_build_args['read_len']) + 'mers.kmer_distrib')

    bracken_build_args_list = [
        '-t', bracken_build_args['threads'],
        '-k', bracken_build_args['kmer_len'],
        '-l', bracken_build_args['read_len'],
        '-d', bracken_build_args['kraken_database'],
    ]

    subprocess.check_call(['bracken-build'] + bracken_build_args_list)

    data_table_entry = {
        "data_tables": {
            data_table_name: [
                {
                    "value": database_value,
                    "name": database_name,
                    "path": database_path,
                }
            ]
        }
    }

    return data_table_entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument('--threads', dest='threads', default=1, help='threads' )
    parser.add_argument('--kmer-len', dest='kmer_len', help='K-mer length' )
    parser.add_argument('--read-len', dest='read_len', help='Read length' )
    parser.add_argument('--kraken-db', dest='kraken_database', help='Kraken Database' )
    parser.add_argument('--database-name', dest='database_name', help='Database Name')
    args = parser.parse_args()

    data_manager_input = json.loads(open(args.data_manager_json).read())

    target_directory = data_manager_input['output_data'][0]['extra_files_path']

    bracken_build_args = {
        'threads': args.threads,
        'kmer_len': args.kmer_len,
        'read_len': args.read_len,
        'kraken_database': args.kraken_database,
    }

    try:
        os.mkdir( target_directory )
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir( target_directory ):
            pass
        else:
            raise

    data_manager_output = {}

    data_manager_output = bracken_build_database(
        target_directory,
        bracken_build_args,
        args.database_name,
    )

    open(args.data_manager_json, 'wb').write(json.dumps(data_manager_output))


if __name__ == "__main__":
    main()
