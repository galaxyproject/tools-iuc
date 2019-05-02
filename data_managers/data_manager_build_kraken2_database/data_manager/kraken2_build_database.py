#!/usr/bin/env python

from __future__ import print_function

import argparse
import datetime
import errno
import json
import os
import shutil
import subprocess
import sys

from enum import Enum

try:
    # Python3
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


DATA_TABLE_NAME = "kraken2_databases"

class KrakenDatabaseTypes(Enum):
    standard = 'standard'
    minikraken = 'minikraken'
    special = 'special'
    custom = 'custom'

    def __str__(self):
        return self.value

class Minikraken2Versions(Enum):
    v1 = 'v1'
    v2 = 'v2'

    def __str__(self):
        return self.value

def kraken2_build_standard(data_manager_dict, kraken2_args, target_directory, data_table_name=DATA_TABLE_NAME):
    now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ")

    database_value = "_".join([
        now,
        "standard",
        "kmer-len", str(kraken2_args["kmer_len"]),
        "minimizer-len", str(kraken2_args["minimizer_len"]),
        "minimizer-spaces", str(kraken2_args["minimizer_spaces"]),
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

    subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--clean',
        '--db', database_path
    ]

    subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    data_table_entry = {
        "value": database_value,
        "name": database_name,
        "path": database_path,
    }

    _add_data_table_entry(data_manager_dict, data_table_entry)


def kraken2_build_minikraken(data_manager_dict, minikraken2_version, target_directory, data_table_name=DATA_TABLE_NAME):

    now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ")

    database_value = "_".join([
        now,
        "minikraken2",
        minikraken2_version,
        "8GB",
    ])

    database_name = " ".join([
        "Minikraken2",
        minikraken2_version,
        "(Created:",
        now + ")"
    ])

    # download the minikraken2 data
    src = urlopen(
        'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_%s_8GB_201904_UPDATE.tgz'
        % minikraken2_version
    )
    with open('tmp_data.tar.gz', 'wb') as dst:
        shutil.copyfileobj(src, dst)
    # unpack the downloaded archive to the target directory
    with tarfile.open('tmp_data.tar.gz', 'r:gz') as fh:
        fh.extractall(target_directory)

    data_table_entry = {
        "value": database_value,
        "name": database_name,
        "path": database_value,
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
    parser.add_argument('--kmer-len', dest='kmer_len', type=int, default=35, help='kmer length')
    parser.add_argument('--minimizer-len', dest='minimizer_len', type=int, default=31, help='minimizer length')
    parser.add_argument('--minimizer-spaces', dest='minimizer_spaces', default=6, help='minimizer spaces')
    parser.add_argument('--threads', dest='threads', default=1, help='threads')
    parser.add_argument('--database-type', dest='database_type', type=KrakenDatabaseTypes, choices=list(KrakenDatabaseTypes), required=True, help='type of kraken database to build')
    parser.add_argument( '--minikraken2-version', dest='minikraken2_version', type=Minikraken2Versions, choices=list(Minikraken2Versions), help='MiniKraken2 version' )
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

    print(args.database_type)
    if str(args.database_type) == 'standard':
        kraken2_args = {
            "kmer_len": args.kmer_len,
            "minimizer_len": args.minimizer_len,
            "minimizer_spaces": args.minimizer_spaces,
            "threads": args.threads,
        }
        kraken2_build_standard(
            data_manager_output,
            kraken2_args,
            target_directory,
        )
    elif str(args.database_type) == 'minikraken':
        kraken2_build_minikraken(
            data_manager_output,
            str(args.minikraken2_version),
            target_directory
        )
    else:
        sys.exit("Invalid database type")

    open(args.data_manager_json, 'w').write(json.dumps(data_manager_output))


if __name__ == "__main__":
    main()
