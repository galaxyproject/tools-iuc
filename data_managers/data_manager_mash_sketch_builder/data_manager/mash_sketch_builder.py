#!/usr/bin/env python

import argparse
import errno
import json
import os
import subprocess
import sys


DATA_TABLE_NAME = "mash_sketches"


def mash_sketch(mash_sketch_args, sketch_name, target_directory, data_table_name=DATA_TABLE_NAME):
    UUID = str(uuid.uuid4())

    os.mkdir(os.path.join(target_directory, UUID))

    sketch_path = os.path.join(UUID, "sketch"),

    args = [
        '--threads', str(kraken2_args["threads"]),
        '-k', str(mash_sketch_args["kmer_size"]),
        '-s', str(mash_sketch_args["sketch_size"]),
        '-o', sketch_path
    ]

    subprocess.check_call(['mash', 'sketch'] + args, cwd=target_directory)

    if kraken2_args["clean"]:
        args = [
            '--threads', str(kraken2_args["threads"]),
            '--clean',
            '--db', database_path
        ]

        subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    data_table_entry = {
        'data_tables': {
            data_table_name: [
                {
                    "value": UUID,
                    "name": sketch_name,
                    "path": sketch_path,
                }
            ]
        }
    }

    return data_table_entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument('--kmer-size', dest='kmer_size', type=int, default=35, help='kmer length')
    parser.add_argument('--sketch-size', dest='sketch_size', type=int, default=31, help='minimizer length')
    parser.add_argument('--threads', dest='threads', default=1, help='threads')
    parser.add_argument('--sketch-name', dest='sketch_name', help='Name for sketch')
    args = parser.parse_args()

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

    mash_sketch_args = {
        "kmer_size": args.kmer_len,
        "sketch_size": args.minimizer_len,
        "threads": args.threads,
    }

    open(args.data_manager_json, 'w').write(json.dumps(data_manager_output, sort_keys=True))


if __name__ == "__main__":
    main()
