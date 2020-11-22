#!/usr/bin/env python

from __future__ import print_function

import argparse
import errno
import json
import os
import subprocess
import uuid


DATA_TABLE_NAME = "kma_index"


def kma_build_index(kma_index_args, index_name, target_directory, data_table_name=DATA_TABLE_NAME):
    UUID = str(uuid.uuid4())

    os.mkdir(os.path.join(target_directory, UUID))

    args = [
        '-k', str(kma_index_args["k"]),
        '-k_t', str(kma_index_args["k_t"]),
        '-k_i', str(kma_index_args["k_i"]),
        '-ML', str(kma_index_args["ML"]),
        '-ht', str(kma_index_args["ht"]),
        '-hq', str(kma_index_args["hq"]),
        '-o', os.path.join(UUID, "index"),
        '-i', " ".join(kma_index_args["fasta"]),
    ]

    subprocess.check_call(' '.join(['kma index'] + args), cwd=target_directory, shell=True)

    data_table_entry = {
        'data_tables': {
            data_table_name: [
                {
                    "value": UUID,
                    "name": index_name,
                    "path": os.path.join(UUID, "index"),
                }
            ]
        }
    }

    return data_table_entry


def main(args):
    with open(args.data_manager_json) as fh:
        data_manager_input = json.load(fh)

    target_directory = data_manager_input['output_data'][0]['extra_files_path']

    try:
        os.mkdir(target_directory)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(target_directory):
            pass
        else:
            raise

    data_manager_output = {}

    kma_index_args = {
        "k": args.k,
        "k_t": args.k_t,
        "k_i": args.k_i,
        "ML": args.ML,
        "ht": args.ht,
        "hq": args.hq,
        "fasta": args.fasta,
    }

    data_manager_output = kma_build_index(
        kma_index_args,
        args.index_name,
        target_directory,
    )

    with open(args.data_manager_json, 'w') as fh:
        json.dump(data_manager_output, fh, sort_keys=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument('--k', dest='k', type=int, default=16, help='')
    parser.add_argument('--k_t', dest='k_t', type=int, default=16, help='')
    parser.add_argument('--k_i', dest='k_i', type=int, default=16, help='')
    parser.add_argument('--ML', dest='ML', type=int, default=16, help='')
    parser.add_argument('--ht', dest='ht', type=float, default=1.0, help='')
    parser.add_argument('--hq', dest='hq', type=float, default=1.0, help='')
    parser.add_argument('--name', dest='index_name', help='')
    parser.add_argument('fasta', nargs='+', help='fasta file(s) to index')
    args = parser.parse_args()
    main(args)
