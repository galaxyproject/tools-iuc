#!/usr/bin/env python

import argparse
import os
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('--db_version', action='store', dest='db_version', help='Version of DRAM databases')
parser.add_argument('--skip_uniref', action='store_true', dest='skip_uniref', default=False, help='Flag to Download and process uniref')
parser.add_argument('--galaxy_data_manager_data_path', action='store', dest='galaxy_data_manager_data_path', help='Absolute Galaxy data manager data path')
parser.add_argument('--output', action='store', dest='output', help='Output file')

args = parser.parse_args()


def get_new_dram_config_entry(db_version, old_entry, new_base_path):
    # Example old_entry:
    # KOfam db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/kofam_profiles.hmm
    base_path, file_name = os.path.split(old_entry)
    # The new entry must be GALAXY_DATA_MANAGER_DATA_PATH/DRAM/${value}/file_name
    return os.path.join(new_base_path, 'DRAM', db_version, file_name)


# At this point the DRAM config will look something like this.
# Processed search databases
# KEGG db: None
# KOfam db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/kofam_profiles.hmm
# KOfam KO list: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/kofam_ko_list.tsv
# UniRef db: None
# Pfam db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/pfam.mmspro
# dbCAN db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/dbCAN-HMMdb-V10.txt
# RefSeq Viral db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/refseq_viral.20220707.mmsdb
# MEROPS peptidase db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/peptidases.20220707.mmsdb
# VOGDB db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/vog_latest_hmms.txt
#
# Descriptions of search database entries
# Pfam hmm dat: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/Pfam-A.hmm.dat.gz
# dbCAN family activities: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/CAZyDB.07292021.fam-activities.txt
# VOG annotations: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/vog_annotations_latest.tsv.gz
#
# Description db: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/description_db.sqlite
#
# DRAM distillation sheets
# Genome summary form: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/genome_summary_form.20220707.tsv
# Module step form: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/module_step_form.20220707.tsv
# ETC module database: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/etc_mdoule_database.20220707.tsv
# Function heatmap form: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/function_heatmap_form.20220707.tsv
# AMG database: /home/galaxies/gvk/jwd/003/3045/working/dataset_4268_files/amg_database.20220707.tsv

# Write the current DRAM CONFIG to a file for processing.
cmd = 'DRAM-setup.py print_config > dram_config.txt'
subprocess.check_call(cmd, shell=True)

# Update the database locations that DRAM sets in it's CONFIG
# to point to the configured GALAXY_DATA_MANAGER_DATA_PATH location
# for the DRAM databases.
cmd = 'DRAM-setup.py set_database_locations'
with open('dram_config.txt', 'r') as fh:
    for line in fh:
        line = line.rstrip('\r\n')
        if line.startswith('KOfam db:'):
            cmd = '%s --kofam_hmm_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('KOfam KO list:'):
            cmd = '%s --kofam_ko_list_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('UniRef db:'):
            if not args.skip_uniref:
                cmd = '%s --uniref_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Pfam db:'):
            cmd = '%s --pfam_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('dbCAN db:'):
            cmd = '%s --dbcan_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('RefSeq Viral db:'):
            cmd = '%s --viral_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('MEROPS peptidase db:'):
            cmd = '%s --peptidase_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('VOGDB db:'):
            cmd = '%s --vogdb_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Pfam hmm dat:'):
            cmd = '%s --pfam_hmm_dat %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('dbCAN family activities:'):
            cmd = '%s --dbcan_fam_activities %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('VOG annotations:'):
            cmd = '%s --vog_annotations %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Description db:'):
            cmd = '%s --description_db_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Genome summary form:'):
            cmd = '%s --genome_summary_form_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Module step form:'):
            cmd = '%s --module_step_form_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('ETC module database:'):
            cmd = '%s --etc_module_database_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('Function heatmap form:'):
            cmd = '%s --function_heatmap_form_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
        elif line.startswith('AMG database:'):
            cmd = '%s --amg_database_loc %s' % (cmd, get_new_dram_config_entry(args.db_version, line, args.galaxy_data_manager_data_path))
cmd = '%s --update_description_db' % cmd
subprocess.check_call(cmd, shell=True)

# Write the new DRAM CONFIG to a file to the output.
cmd = 'DRAM-setup.py print_config > %s' % args.output
subprocess.check_call(cmd, shell=True)
