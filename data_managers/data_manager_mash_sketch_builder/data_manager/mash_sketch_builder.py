#!/usr/bin/env python

import argparse
import errno
import json
import os
import subprocess
import uuid


DATA_TABLE_NAME = "mash_sketches"


def mash_sketch(mash_sketch_args, sketch_name, target_directory, data_table_name=DATA_TABLE_NAME):
    UUID = str(uuid.uuid4())

    os.mkdir(os.path.join(target_directory, UUID))

    sketch_path = os.path.join(target_directory, UUID, "sketch")

    args = [
        '-k', str(mash_sketch_args["kmer_size"]),
        '-s', str(mash_sketch_args["sketch_size"]),
        '-w', str(mash_sketch_args["probability_threshold"]),
        '-o', str(sketch_path),
        '-p', str(mash_sketch_args["threads"]),
        str(mash_sketch_args["fasta"]),
    ]

    if mash_sketch_args["individual_sequences"]:
        args = args + ["-i"]

    subprocess.check_call(['mash', 'sketch'] + args, cwd=target_directory)

    data_table_entry = {
        'data_tables': {
            data_table_name: [
                {
                    "value": UUID,
                    "name": sketch_name,
                    "path": UUID,
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
    parser.add_argument('--probability-threshold', dest='probability_threshold', type=float, default=0.01, help='Probability threshold for warning about low k-mer size')
    parser.add_argument('--individual-sequences', dest='individual_sequences', action='store_true', default=False, help='Sketch individual sequences (for multi-fasta files)')
    parser.add_argument('--fasta', dest='fasta', help='Fasta file to sketch')
    parser.add_argument('--threads', dest='threads', default=1, help='threads')
    parser.add_argument('--sketch-name', dest='sketch_name', help='Name for sketch')
    args = parser.parse_args()

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

    mash_sketch_args = {
        "kmer_size": args.kmer_size,
        "sketch_size": args.sketch_size,
        "probability_threshold": args.probability_threshold,
        "fasta": args.fasta,
        "individual_sequences": args.individual_sequences,
        "threads": args.threads,
    }

    data_manager_output = mash_sketch(
        mash_sketch_args,
        args.sketch_name,
        target_directory,
    )

    with open(args.data_manager_json, 'w') as fh:
        json.dump(data_manager_output, fh, sort_keys=True)


if __name__ == "__main__":
    main()
