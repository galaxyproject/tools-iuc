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
import tarfile
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


class SpecialDatabaseTypes(Enum):
    rdp = 'rdp'
    greengenes = 'greengenes'
    silva = 'silva'

    def __str__(self):
        return self.value


class Minikraken2Versions(Enum):
    v1 = 'v1'
    v2 = 'v2'

    def __str__(self):
        return self.value


def kraken2_build_standard(kraken2_args, target_directory, data_table_name=DATA_TABLE_NAME):
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
        'data_tables': {
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


def kraken2_build_minikraken(minikraken2_version, target_directory, data_table_name=DATA_TABLE_NAME):

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

    database_path = database_value

    # download the minikraken2 data
    src = urlopen(
        'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_%s_8GB_201904_UPDATE.tgz'
        % minikraken2_version
    )
    with open('tmp_data.tar.gz', 'wb') as dst:
        shutil.copyfileobj(src, dst)
    # unpack the downloaded archive to the target directory
    with tarfile.open('tmp_data.tar.gz', 'r:gz') as fh:
        for member in fh.getmembers():
            if member.isreg():
                member.name = os.path.basename(member.name)
                fh.extract(member, os.path.join(target_directory, database_path))


    data_table_entry = {
        'data_tables': {
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


def kraken2_build_special(kraken2_args, target_directory, data_table_name=DATA_TABLE_NAME):

    now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H%M%SZ")

    special_database_names = {
        "rdp": "RDP",
        "greengenes": "Greengenes",
        "silva": "Silva",
    }

    database_value = "_".join([
        now,
        kraken2_args["special_database_type"],
        "kmer-len", str(kraken2_args["kmer_len"]),
        "minimizer-len", str(kraken2_args["minimizer_len"]),
        "minimizer-spaces", str(kraken2_args["minimizer_spaces"]),
    ])

    database_name = " ".join([
        special_database_names[kraken2_args["special_database_type"]],
        "(Created:",
        now + ",",
        "kmer-len=" + str(kraken2_args["kmer_len"]) + ",",
        "minimizer-len=" + str(kraken2_args["minimizer_len"]) + ",",
        "minimizer-spaces=" + str(kraken2_args["minimizer_spaces"]) + ")",
    ])

    database_path = database_value

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--special', kraken2_args["special_database_type"],
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
        'data_tables': {
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


def kraken2_build_custom(kraken2_args, custom_database_name, target_directory, data_table_name=DATA_TABLE_NAME):

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--download-taxonomy',
        '--db', custom_database_name
    ]

    subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--add-to-library', kraken2_args["custom_fasta"],
        '--db', custom_database_name
    ]

    subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--build',
        '--kmer-len', str(kraken2_args["kmer_len"]),
        '--minimizer-len', str(kraken2_args["minimizer_len"]),
        '--minimizer-spaces', str(kraken2_args["minimizer_spaces"]),
        '--db', custom_database_name
    ]

    subprocess.check_call(['kraken2-build'] + args, cwd=target_directory)

    args = [
        '--threads', str(kraken2_args["threads"]),
        '--clean',
        '--db', custom_database_name
    ]

    subprocess.check_call(['kraken2-build'] + args, target_directory)

    data_table_entry = {
        'data_tables': {
            data_table_name: [
                {
                    "value": custom_database_name,
                    "name": custom_database_name,
                    "path": custom_database_name
                }
            ]
        }
    }

    return data_table_entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('data_manager_json')
    parser.add_argument('--kmer-len', dest='kmer_len', type=int, default=35, help='kmer length')
    parser.add_argument('--minimizer-len', dest='minimizer_len', type=int, default=31, help='minimizer length')
    parser.add_argument('--minimizer-spaces', dest='minimizer_spaces', default=6, help='minimizer spaces')
    parser.add_argument('--threads', dest='threads', default=1, help='threads')
    parser.add_argument('--database-type', dest='database_type', type=KrakenDatabaseTypes, choices=list(KrakenDatabaseTypes), required=True, help='type of kraken database to build')
    parser.add_argument('--minikraken2-version', dest='minikraken2_version', type=Minikraken2Versions, choices=list(Minikraken2Versions), help='MiniKraken2 version (only applies to --database-type minikraken)')
    parser.add_argument('--special-database-type', dest='special_database_type', type=SpecialDatabaseTypes, choices=list(SpecialDatabaseTypes), help='type of special database to build (only applies to --database-type special)')
    parser.add_argument('--custom-fasta', dest='custom_fasta', help='fasta file for custom database (only applies to --database-type custom)')
    parser.add_argument('--custom-database-name', dest='custom_database_name', help='Name for custom database (only applies to --database-type custom)')
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

    if str(args.database_type) == 'standard':
        kraken2_args = {
            "kmer_len": args.kmer_len,
            "minimizer_len": args.minimizer_len,
            "minimizer_spaces": args.minimizer_spaces,
            "threads": args.threads,
        }
        data_manager_output = kraken2_build_standard(
            kraken2_args,
            target_directory,
        )
    elif str(args.database_type) == 'minikraken':
        data_manager_output = kraken2_build_minikraken(
            str(args.minikraken2_version),
            target_directory
        )
    elif str(args.database_type) == 'special':
        kraken2_args = {
            "special_database_type": str(args.special_database_type),
            "kmer_len": args.kmer_len,
            "minimizer_len": args.minimizer_len,
            "minimizer_spaces": args.minimizer_spaces,
            "threads": args.threads,
        }
        data_manager_output = kraken2_build_special(
            kraken2_args,
            target_directory,
        )
    elif str(args.database_type) == 'custom':
        kraken2_args = {
            "custom_fasta": args.custom_fasta,
            "kmer_len": args.kmer_len,
            "minimizer_len": args.minimizer_len,
            "minimizer_spaces": args.minimizer_spaces,
            "threads": args.threads,
        }
        data_manager_output = kraken2_build_custom(
            kraken2_args,
            args.custom_database_name,
            target_directory,
        )
    else:
        sys.exit("Invalid database type")

    open(args.data_manager_json, 'w').write(json.dumps(data_manager_output))


if __name__ == "__main__":
    main()
