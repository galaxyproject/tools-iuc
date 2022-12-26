#!/usr/bin/env python

import json
import optparse
import os
import subprocess
import sys


def get_id_name(params, dbkey, fasta_description=None):
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey

    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = fasta_description
        if not sequence_name:
            sequence_name = dbkey
    return sequence_id, sequence_name


def build_malt_index(data_manager_dict, fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, sequence_type, shapes, max_hits_per_seed, protein_reduct):
    # The malt-build program produces a directory of files,
    # so the data table path entry will be a directory and
    # not an index file.
    fasta_base_name = os.path.split(fasta_filename)[-1]
    sym_linked_fasta_filename = os.path.join(target_directory, fasta_base_name)
    os.symlink(fasta_filename, sym_linked_fasta_filename)
    args = ['malt-build', '--input', sym_linked_fasta_filename, '--sequenceType', sequence_type, '--index', target_directory]
    threads = os.environ.get('GALAXY_SLOTS')
    if threads:
        args.extend(['--threads', threads])
    if shapes is not None:
        args.extend(['--shapes', shapes])
    if max_hits_per_seed is not None:
        args.extend(['--maxHitsPerSeed', max_hits_per_seed])
    if protein_reduct is not None:
        args.extend(['--proteinReduct', protein_reduct])
    proc = subprocess.Popen(args=args, shell=False, cwd=target_directory)
    return_code = proc.wait()
    if return_code:
        sys.exit('Error building index, return_code: %d' % return_code)
    # Remove unwanted files from the output directory.
    os.remove(sym_linked_fasta_filename)
    # The path entry here is the directory
    # where the index files will be located,
    # not a single index file (malt-build
    # produces a directory if files, which
    # is considered an index..
    data_table_entry = dict(value=sequence_id, dbkey=dbkey, name=sequence_name, path=None)
    _add_data_table_entry(data_manager_dict, data_table_entry)


def _add_data_table_entry(data_manager_dict, data_table_entry):
    data_table_name = "malt_indices"
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get(data_table_name, [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def main():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--fasta_filename', dest='fasta_filename', action='store', type="string", help='fasta filename')
    parser.add_option('-d', '--fasta_dbkey', dest='fasta_dbkey', action='store', type="string", help='fasta dbkey')
    parser.add_option('-t', '--fasta_description', dest='fasta_description', action='store', type="string", default=None, help='fasta description')
    parser.add_option('-e', '--sequence_type', dest='sequence_type', action='store', type="string", help='DNA or Protein sequences')
    parser.add_option('-p', '--shapes', dest='shapes', action='store', type="string", default=None, help='Comma-separated list of seed shapes')
    parser.add_option('-m', '--max_hits_per_seed', dest='max_hits_per_seed', action='store', type="string", default=None, help='Maximum number of hits per seed')
    parser.add_option('-r', '--protein_reduct', dest='protein_reduct', action='store', type="string", default=None, help='Name or definition of protein alphabet reduction')
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

    # Build the index.
    build_malt_index(data_manager_dict, options.fasta_filename, params, target_directory, dbkey, sequence_id, sequence_name, options.sequence_type, options.shapes, options.max_hits_per_seed, options.protein_reduct)

    # Save info to json file.
    with open(filename, 'w') as fh:
        json.dump(data_manager_dict, fh, sort_keys=True)


if __name__ == "__main__":
    main()
