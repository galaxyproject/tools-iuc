#!/usr/bin/env python

import argparse
import json
import os
import subprocess


parser = argparse.ArgumentParser()

parser.add_argument('--kofam_hmm_loc', action='store', dest='kofam_hmm_loc', default=None, help='hmm file for KOfam')
parser.add_argument('--kofam_ko_list_loc', action='store', dest='kofam_ko_list_loc', default=None, help='KOfam ko list file')
parser.add_argument('--skip_uniref', action='store', dest='skip_uniref', default='no', help='Flag to Download and process uniref')
parser.add_argument('--uniref_loc', action='store', dest='uniref_loc', default=None, help='uniref file')
parser.add_argument('--uniref_version', action='store', dest='uniref_version', type=int, default=90, help='uniref version to download')
parser.add_argument('--pfam_loc', action='store', dest='pfam_loc', default=None, help='pfam-A full file')
parser.add_argument('--pfam_hmm_dat', action='store', dest='pfam_hmm_dat', help='pfam hmm .dat file to get PF descriptions')
parser.add_argument('--dbcan_loc', action='store', dest='dbcan_loc', default=None, help='dbCAN file')
parser.add_argument('--dbcan_fam_activities', action='store', dest='dbcan_fam_activities', default=None, help='CAZY family activities file')
parser.add_argument('--dbcan_version', action='store', dest='dbcan_version', type=int, default=10, help='Version of dbCAN to use')
parser.add_argument('--vogdb_loc', action='store', dest='vogdb_loc', default=None, help='hmm file for vogdb')
parser.add_argument('--vog_annotations', action='store', dest='vog_annotations', default=None, help='vogdb annotations file')
parser.add_argument('--viral_loc', action='store', dest='viral_loc', default=None, help='merged viral protein faa file')
parser.add_argument('--peptidase_loc', action='store', dest='peptidase_loc', default=None, help='MEROPS peptidase fasta file')
parser.add_argument('--genome_summary_form_loc', action='store', dest='genome_summary_form_loc', default=None, help='genome summary form file')
parser.add_argument('--module_step_form_loc', action='store', dest='module_step_form_loc', default=None, help='module step form file')
parser.add_argument('--etc_module_database_loc', action='store', dest='etc_module_database_loc', default=None, help='etc module database file')
parser.add_argument('--function_heatmap_form_loc', action='store', dest='function_heatmap_form_loc', default=None, help='function heatmap form file')
parser.add_argument('--amg_database_loc', action='store', dest='amg_database_loc', default=None, help='amg database file')
parser.add_argument('--db_version', action='store', dest='db_version', help='Version of DRAM databases')
parser.add_argument('--threads', action='store', dest='threads', type=int, help='Number of processes')
parser.add_argument('--out_file', action='store', dest='out_file', help='JSON output file')

args = parser.parse_args()

with open(args.out_file) as fh:
    params = json.load(fh)

target_directory = params['output_data'][0]['extra_files_path']
os.makedirs(target_directory)

# Download the data.
cmd = 'DRAM-setup.py prepare_databases --output_dir %s' % target_directory
if args.kofam_hmm_loc is not None:
    cmd = '%s --kofam_hmm_loc %s' % (cmd, args.kofam_hmm_loc)
if args.kofam_ko_list_loc is not None:
    cmd = '%s --kofam_ko_list_loc %s' % (cmd, args.kofam_ko_list_loc)
if args.skip_uniref == 'yes':
    cmd = '%s --skip_uniref' % cmd
else:
    if args.uniref_loc is not None:
        cmd = '%s --uniref_loc %s' % (cmd, args.uniref_loc)
    cmd = '%s --uniref_version %d' % (cmd, args.uniref_version)
if args.pfam_loc is not None:
    cmd = '%s --pfam_loc %s' % (cmd, args.pfam_loc)
if args.pfam_hmm_dat is not None:
    cmd = '%s --pfam_hmm_dat %s' % (cmd, args.pfam_hmm_dat)
if args.dbcan_loc is not None:
    cmd = '%s --dbcan_loc %s' % (cmd, args.dbcan_loc)
if args.dbcan_fam_activities is not None:
    cmd = '%s --dbcan_fam_activities %s' % (cmd, args.dbcan_fam_activities)
cmd = '%s --dbcan_version %d' % (cmd, args.dbcan_version)
if args.vogdb_loc is not None:
    cmd = '%s --vogdb_loc %s' % (cmd, args.vogdb_loc)
if args.vog_annotations is not None:
    cmd = '%s --vog_annotations %s' % (cmd, args.vog_annotations)
if args.viral_loc is not None:
    cmd = '%s --viral_loc %s' % (cmd, args.viral_loc)
if args.peptidase_loc is not None:
    cmd = '%s --peptidase_loc %s' % (cmd, args.peptidase_loc)
if args.genome_summary_form_loc is not None:
    cmd = '%s --genome_summary_form_loc %s' % (cmd, args.genome_summary_form_loc)
if args.module_step_form_loc is not None:
    cmd = '%s --module_step_form_loc %s' % (cmd, args.module_step_form_loc)
if args.etc_module_database_loc is not None:
    cmd = '%s --etc_module_database_loc %s' % (cmd, args.etc_module_database_loc)
if args.function_heatmap_form_loc is not None:
    cmd = '%s --function_heatmap_form_loc %s' % (cmd, args.function_heatmap_form_loc)
if args.amg_database_loc is not None:
    cmd = '%s --amg_database_loc %s' % (cmd, args.amg_database_loc)
cmd = '%s --threads %d' % (cmd, args.threads)

subprocess.check_call(cmd, shell=True)

data_manager_json = {'data_tables': {}}
data_manager_entry = {}
data_manager_entry['value'] = args.db_version
data_manager_entry['name'] = 'DRAM %s databases' % args.db_version
data_manager_entry['path'] = target_directory
data_manager_json['data_tables']['dram_databases'] = data_manager_entry

with open(args.out_file, 'w') as fh:
    json.dump(data_manager_json, fh, sort_keys=True)
