#!/usr/bin/env python

from __future__ import print_function

import argparse
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


def kraken2_build(data_manager_dict, kraken2_args, database_name, params, target_directory, data_table_name=DATA_TABLE_NAME):

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--download-taxonomy',
        '--db', database_name
    ]

    run(['kraken2-build'] + args, target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--add-to-library', kraken2_args["fasta"],
        '--db', database_name
    ]

    run(['kraken2-build'] + args, target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--build',
        '--kmer-len', str(kraken2_args["kmer_len"]),
        '--minimizer-len', str(kraken2_args["minimizer_len"]),
        '--minimizer-spaces', str(kraken2_args["minimizer_spaces"]),
        '--db', database_name
    ]

    run(['kraken2-build'] + args, target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--clean',
        '--db', database_name
    ]

    run(['kraken2-build'] + args, target_directory)

    data_table_entry = {
        "value": database_name,
        "name": database_name,
        "path": database_name
    }

    _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry)


def _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('params')
    parser.add_argument( '-d', '--db', dest='database_name', help='database name' )
    parser.add_argument( '-k', '--kmer-len', dest='kmer_len', type=int, default=35, help='kmer length' )
    parser.add_argument( '-m', '--minimizer-len', dest='minimizer_len', type=int, default=31, help='minimizer length' )
    parser.add_argument( '-s', '--minimizer-spaces', dest='minimizer_spaces', default=6, help='minimizer spaces' )
    parser.add_argument( '-f', '--fasta', dest='fasta', help='fasta' )
    parser.add_argument( '-t', '--threads', dest='threads', default=1, help='threads' )
    args = parser.parse_args()

    kraken2_args = {
        "kmer_len": args.kmer_len,
        "minimizer_len": args.minimizer_len,
        "minimizer_spaces": args.minimizer_spaces,
        "fasta": args.fasta,
        "threads": args.threads,
    }

    params = json.loads(open(args.params).read())
    target_directory = params['output_data'][0]['extra_files_path']

    try:
        os.mkdir( target_directory )
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir( target_directory ):
            pass
        else:
            raise

    data_manager_dict = {}

    # build the index
    kraken2_build(
        data_manager_dict,
        kraken2_args,
        params,
        target_directory
    )

    # save info to json file
    open(args.params, 'wb').write(json.dumps(data_manager_dict))


if __name__ == "__main__":
    main()
