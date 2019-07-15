#!/usr/bin/env python
# Dan Blankenberg
from __future__ import print_function

import argparse
import os
import subprocess
import sys

import yaml

try:
    from shlex import quote as cmd_quote
except ImportError:
    from pipes import quote as cmd_quote


def build_bowtie2_index(fasta_filename, target_directory, index_id):
    # TODO: allow multiple FASTA input files
    fasta_base_name = os.path.split( fasta_filename )[-1]
    sym_linked_fasta_filename = os.path.join( target_directory, fasta_base_name )
    os.symlink( fasta_filename, sym_linked_fasta_filename )
    args = ['bowtie2-build', sym_linked_fasta_filename, index_id]
    proc = subprocess.Popen(args=args, shell=False, cwd=target_directory)
    return_code = proc.wait()
    if return_code:
        print("Error building index.", file=sys.stderr)
        sys.exit(return_code)
    return [' '.join(cmd_quote(arg) for arg in args)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata')
    parser.add_argument( '-t', '--target-directory', default=None)
    args = parser.parse_args()

    params = yaml.safe_load(open(args.metadata))
    params['cmds'] = []
    if args.target_directory:
        if not os.path.isdir(args.target_directory):
            os.mkdir(args.target_directory)
    else:
        args.target_directory = os.getcwd()

    dbkey = params['dbkey']
    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )

    # build the index
    params['cmds'] += build_bowtie2_index(params['source'], args.target_directory, params['id'])
    with open(args.metadata, 'w') as fo:
        yaml.dump(params, fo, allow_unicode=False, default_flow_style=False)


if __name__ == "__main__":
    main()
