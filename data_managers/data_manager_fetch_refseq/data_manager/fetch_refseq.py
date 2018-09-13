#!/usr/bin/env python

from __future__ import division, print_function

import argparse
import functools
import gzip
import json
import os
import os.path
import sys
from datetime import date
from multiprocessing import Process, Queue

import requests

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO
# Refseq structure
# - Release number
# - Divisions
#   1. archea
#   2. bacteria
#   3. fungi
#   4. invertebrate
#   5. mitochondrion
#   6. other
#   7. plant
#   8. plasmid
#   9. plastid
#  10. protozoa
#  11. vertebrate mammalian
#  12. vertebrate other
#  13. viral
# within each division
# DIVNAME.\d+(.\d+)?.(genomic|protein|rna).(fna|gbff|faa|gpff).gz
#  where fna and faa are FASTA, gbff and gpff are Genbank


def _add_data_table_entry(data_manager_dict, data_table_entry, data_table_name):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get('all_fasta', [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def unzip_to(conn, out_dir, output_filename, chunk_size=4096, debug=False, compress=False):
    input_filename = conn.get()
    if compress:
        open_output = gzip.open
    else:
        open_output = open
    with open_output(os.path.join(out_dir, output_filename), 'wb') as output_file:
        while input_filename != 'STOP':
            if debug:
                print('Reading', input_filename, file=sys.stderr)
            with gzip.open(input_filename, 'rb') as input_file:
                read_chunk = functools.partial(input_file.read, (chunk_size))
                for data in iter(read_chunk, b''):  # use b'' as a sentinel to stop the loop. note '' != b'' in Python 3
                    output_file.write(data)
            os.unlink(input_filename)
            input_filename = conn.get()


def get_refseq_division(division_name, mol_types, output_directory, debug=False, compress=False):
    base_url = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/'
    valid_divisions = set(['archea', 'bacteria', 'complete', 'fungi', 'invertebrate', 'mitochondrion', 'other',
                          'plant', 'plasmid', 'plastid', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral'])
    ending_mappings = {
        'genomic': '.genomic.fna.gz',
        'protein': '.protein.faa.gz',
        'rna': 'rna.fna.gz'
    }
    assert division_name in valid_divisions, "Unknown division name ({})".format(division_name)
    for mol_type in mol_types:
        assert mol_type in ending_mappings, "Unknown molecule type ({})".format(mol_type)
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    release_num_file = base_url + 'RELEASE_NUMBER'
    r = requests.get(release_num_file)
    release_num = r.text.strip()
    division_base_url = base_url + division_name
    if debug:
        print('Retrieving {}'.format(division_base_url), file=sys.stderr)
    r = requests.get(division_base_url)
    listing_text = r.text

    unzip_queues = {}
    unzip_processes = []
    final_output_filenames = []
    for mol_type in mol_types:
        q = unzip_queues[mol_type] = Queue()
        output_filename = division_name + '.' + release_num + '.' + mol_type + '.fasta'
        if compress:
            output_filename += '.gz'
        final_output_filenames.append(output_filename)
        unzip_processes.append(Process(target=unzip_to, args=(q, output_directory, output_filename),
                                       kwargs=dict(debug=debug, compress=compress)))
        unzip_processes[-1].start()

    # sample line: <a href="vertebrate_other.86.genomic.gbff.gz">vertebrate_other.86.genomic.gbff.gz</a>   2018-07-13 00:59   10M
    for line in StringIO(listing_text):
        if '.gz' not in line:
            continue
        parts = line.split('"')
        assert len(parts) == 3, "Unexpected line format: {}".format(line.rstrip())
        filename = parts[1]
        for mol_type in mol_types:
            ending = ending_mappings[mol_type]
            if filename.endswith(ending):
                if debug:
                    print('Downloading:', filename, ending, mol_type, file=sys.stderr)
                output_filename = os.path.join(output_directory, filename)
                with open(output_filename, 'wb') as output_file:
                    r = requests.get(division_base_url + '/' + filename)
                    for chunk in r.iter_content(chunk_size=4096):
                        output_file.write(chunk)
                conn = unzip_queues[mol_type]
                conn.put(output_filename)

    for mol_type in mol_types:
        conn = unzip_queues[mol_type]
        conn.put('STOP')

    return [release_num, final_output_filenames]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download RefSeq databases')
    parser.add_argument('--debug', default=False, action='store_true', help='Print debugging output to stderr (verbose)')
    parser.add_argument('--compress', default=False, action='store_true', help='Compress output files')
    parser.add_argument('--output_directory', default='tmp', help='Directory to write output to')
    parser.add_argument('--galaxy_datamanager_filename', help='Galaxy JSON format file describing data manager inputs')
    parser.add_argument('--division_names', help='RefSeq divisions to download')
    parser.add_argument('--mol_types', help='Molecule types (genomic, rna, protein) to fetch')
    parser.add_argument('--pin_date', help='Force download date to this version string')
    args = parser.parse_args()

    division_names = args.division_names.split(',')
    mol_types = args.mol_types.split(',')
    if args.galaxy_datamanager_filename is not None:
        dm_opts = json.loads(open(args.galaxy_datamanager_filename).read())
        output_directory = dm_opts['output_data'][0]['extra_files_path']  # take the extra_files_path of the first output parameter
        data_manager_dict = {}
    else:
        output_directory = args.output_directory
    for division_name in division_names:
        if args.pin_date is not None:
            today_str = args.pin_date
        else:
            today_str = date.today().strftime('%Y-%m-%d')  # ISO 8601 date format
        [release_num, fasta_files] = get_refseq_division(division_name, mol_types, output_directory, args.debug, args.compress)
        if args.galaxy_datamanager_filename is not None:
            for i, mol_type in enumerate(mol_types):
                assert mol_type in fasta_files[i], "Filename does not contain expected mol_type ({}, {})".format(mol_type, fasta_files[i])
                unique_key = division_name + '.' + release_num + '.' + mol_type  # note: this is now same as dbkey
                dbkey = division_name + '.' + release_num + '.' + mol_type
                desc = 'RefSeq ' + division_name + ' Release ' + release_num + ' ' + mol_type + ' (' + today_str + ')'
                path = os.path.join(output_directory, fasta_files[i])
                _add_data_table_entry(data_manager_dict=data_manager_dict,
                                      data_table_entry=dict(value=unique_key, dbkey=dbkey, name=desc, path=path),
                                      data_table_name='all_fasta')
            open(args.galaxy_datamanager_filename, 'w').write(json.dumps(data_manager_dict, sort_keys=True))
