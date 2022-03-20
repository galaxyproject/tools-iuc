#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import bz2
import cPickle as pickle
import json


def transform_json_to_pkl(args):
    with open(args.json_input, 'r') as json_file:
        json_str = json_file.read()
        metadata = json.loads(json_str)

        for marker in metadata["markers"]:
            a_set = set(metadata["markers"][marker]["ext"])
            metadata["markers"][marker]["ext"] = a_set

    pkl_output = bz2.BZ2File(args.pkl_output, 'w')
    pickle.dump(metadata, pkl_output, pickle.HIGHEST_PROTOCOL)
    pkl_output.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--json_input', required=True)
    parser.add_argument('--pkl_output', required=True)
    args = parser.parse_args()

    transform_json_to_pkl(args)
