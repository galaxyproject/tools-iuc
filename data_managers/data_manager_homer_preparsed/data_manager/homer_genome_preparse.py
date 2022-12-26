#!/usr/bin/env python
# Dan Blankenberg for bowtie2
# Modified by Lucille Delisle for homer
from __future__ import print_function

import json
import optparse
import os
import subprocess
import sys

DEFAULT_DATA_TABLE_NAME = "homer_preparse"


def get_id_name(params, dbkey, fasta_description=None):
    # TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey

    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = fasta_description
        if not sequence_name:
            sequence_name = dbkey
    return sequence_id, sequence_name


def homer_preparse(data_manager_dict, fasta_filename, params, target_directory, dbkey, sequence_id,
                   sequence_name, size, mask, version,
                   data_table_name=DEFAULT_DATA_TABLE_NAME):
    args = ['preparseGenome.pl', fasta_filename, '-size', str(size), '-preparsedDir', target_directory]
    if mask:
        args.append('-mask')
    proc = subprocess.Popen(args=args, shell=False, cwd=target_directory)
    return_code = proc.wait()
    if return_code:
        print("Error preparsing genome.", file=sys.stderr)
        sys.exit(return_code)
    mask_suffix = 'r' if mask else ''
    mask_suffix_name = ' masked' if mask else ''
    data_table_entry = dict(value=sequence_id + mask_suffix + '_' + str(size), dbkey=dbkey,
                            mask=str(mask), size=str(size), name=sequence_name + mask_suffix_name + ' (' + str(size) + 'bp)',
                            path=sequence_id + mask_suffix + '_' + str(size),
                            path_fasta=fasta_filename,
                            version=version)
    _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry)


def _add_data_table_entry(data_manager_dict, data_table_name, data_table_entry):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get(data_table_name, [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--fasta_filename', dest='fasta_filename', action='store', type="string", default=None, help='fasta_filename')
    parser.add_option('-d', '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", default=None, help='fasta_dbkey')
    parser.add_option('-t', '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta_description')
    parser.add_option('-s', '--size', dest='size', action='store', type="int", default=200, help='fragment size')
    parser.add_option('-m', '--mask', dest='mask', action='store_true', default=False, help='mask the lower case bases (repeats)')
    parser.add_option('-n', '--data_table_name', dest='data_table_name', action='store', type="string", default=None, help='data_table_name')
    parser.add_option('--index_version', dest='index_version', action='store', type="string", default=None, help='index version')
    (options, args) = parser.parse_args()

    filename = args[0]

    with open(filename) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    data_manager_dict = {}

    dbkey = options.fasta_dbkey

    if dbkey in [None, '', '?']:
        raise Exception('"%s" is not a valid dbkey. You must specify a valid dbkey.' % (dbkey))

    sequence_id, sequence_name = get_id_name(params, dbkey=dbkey, fasta_description=options.fasta_description)

    # preparse the genome
    homer_preparse(data_manager_dict, options.fasta_filename, params, target_directory, dbkey, sequence_id,
                   sequence_name, options.size, options.mask, options.index_version,
                   data_table_name=options.data_table_name or DEFAULT_DATA_TABLE_NAME)

    # save info to json file
    with open(filename, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
