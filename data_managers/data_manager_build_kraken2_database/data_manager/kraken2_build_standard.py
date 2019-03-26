#!/usr/bin/env python

from __future__ import print_function

import argparse
import datetime
import errno
import json
import os
import subprocess
import sys


DATA_TABLE_NAME = "kraken2_databases"


def run(args, cwd):
    proc = subprocess.Popen(args=args, shell=False, cwd=cwd)
    return_code = proc.wait()
    if return_code:
        print("Error building database.", file=sys.stderr)
        sys.exit( return_code )

def kraken2_build_standard(data_manager_dict, kraken2_args, target_directory, data_table_name=DATA_TABLE_NAME):
    now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ")

    database_value = "_".join([
        now,
        "standard",
        "kmer-len" + str(kraken2_args["kmer_len"]),
        "minimizer-len" + str(kraken2_args["minimizer_len"]),
        "minimizer-spaces" + str(kraken2_args["minimizer_spaces"]),
    ])

    database_name = " ".join([
        "Standard",
        "(Created:",
        now + ",",
        "kmer-len=" + str(kraken2_args["kmer_len"]) + ",",
        "minimizer-len=" + str(kraken2_args["minimizer_len"]) + ",",
        "minimizer-spaces=" + str(kraken2_args["minimizer_spaces"]) + ")",
    ])

    database_path = database_value

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--standard',
        '--kmer-len', str(kraken2_args["kmer_len"]),
        '--minimizer-len', str(kraken2_args["minimizer_len"]),
        '--minimizer-spaces', str(kraken2_args["minimizer_spaces"]),
        '--db', database_path
    ]

    run(['kraken2-build'] + args, target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--clean',
        '--db', database_path
    ]

    run(['kraken2-build'] + args, target_directory)

    data_table_entry = {
        "value": database_value,
        "name": database_name,
        "path": database_path,
    }

    _add_data_table_entry(data_manager_dict, data_table_entry)


def _add_data_table_entry(data_manager_dict, data_table_entry, data_table_name=DATA_TABLE_NAME):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][data_table_name].append( data_table_entry )
    return data_manager_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument( '-k', '--kmer-len', dest='kmer_len', type=int, default=35, help='kmer length' )
    parser.add_argument( '-m', '--minimizer-len', dest='minimizer_len', type=int, default=31, help='minimizer length' )
    parser.add_argument( '-s', '--minimizer-spaces', dest='minimizer_spaces', default=6, help='minimizer spaces' )
    parser.add_argument( '-t', '--threads', dest='threads', default=1, help='threads' )
    args = parser.parse_args()

    kraken2_args = {
        "kmer_len": args.kmer_len,
        "minimizer_len": args.minimizer_len,
        "minimizer_spaces": args.minimizer_spaces,
        "threads": args.threads,
    }

    data_manager_input = json.loads(open(args.data_manager_json).read())

    target_directory = data_manager_input['output_data'][0]['extra_files_path']

    try:
        os.mkdir( target_directory )
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir( target_directory ):
            pass
        else:
            raise

    data_manager_output = {}

    kraken2_build_standard(
        data_manager_output,
        kraken2_args,
        target_directory,
    )

    open(args.data_manager_json, 'wb').write(json.dumps(data_manager_output))


if __name__ == "__main__":
    main()
