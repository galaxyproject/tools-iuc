#!/usr/bin/env python

import argparse
import json
import os
import sys
import subprocess
import dataclasses


# TypeError: non-default argument 'manifest' follows default argument

@dataclasses.dataclass
class LocFileEntry:
    tag: str
    description: str
    manifest: str
    local_cache_dir: str
    tool_cache_dir: str
    default: int = 0
    phone_home: int = 0

class LocFile:
    def __init__(self, pathname):
        self.pathname = pathname
        self.entries = {}
        self.entries_to_add = {}
        self.entries_to_remove = {}
        self._load_file()

    def _load_file(self):
        with open(self.pathname) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    continue

                tag, description, default, phone_home, manifest, local_cache_dir, tool_cache_dir = line.split('\t', 6)

                try:
                    default = int(default)
                except:
                    default = 0

                try:
                    phone_home = int(phone_home)
                except:
                    phone_home = 0

                entry = LocFileEntry(tag, description, manifest, local_cache_dir, tool_cache_dir, default, phone_home)
                self.entries[entry.tag] = entry

    def add_entry(self, entry):
        if entry.tag in self.entries:
            sys.exit(f'entry with tag {args.tag} already exists')
        else:
            self.entries[entry.tag] = entry
            self.entries_to_add[entry.tag] = entry

    def remove_entry(self, entry):
        if entry.tag not in self.entries:
            sys.exit(f'entry with tag {args.tag} does not exist')
        else:
            self.entries_to_remove[entry.tag] = entry
            del self.entries[entry.tag]

    def get_entry(self, tag):
        return self.entries.get(tag, None)

    def configure_database(self, args, params):
        entry = self.get_entry(args.tag)
        new_entry = LocFileEntry(args.tag, args.description, args.manifest, args.local_cache_dir, args.tool_cache_dir, args.default, args.phone_home)

        if args.mode in ['update', 'remove']:
            self.remove_entry(entry)

        if args.mode in ['add', 'update']:
            self.add_entry(new_entry)

        update_dict = self.get_update_dict()
        with open(args.output_file, 'w') as fh:
            json.dump(update_dict, fh, indent=2, sort_keys=True)

    def get_update_dict(self):
        update_dict = {
            'data_tables': {
                'ncbi_fcs_gx_databases': {
                }
            }
        }

        database_dict = update_dict['data_tables']['ncbi_fcs_gx_databases']

        for tag in self.entries_to_add.keys():
            entry = self.entries_to_add[tag]
            database_dict.setdefault('add', [])

            entry_dict = dataclasses.asdict(entry)

            entry_dict['value'] = entry_dict['tag']
            del entry_dict['tag']

            entry_dict['name'] = entry_dict['description']
            del entry_dict['description']

            if 'default' in entry_dict:
                entry_dict['default'] = str(entry_dict['default'])
            if 'phone_home' in entry_dict:
                entry_dict['phone_home'] = str(entry_dict['phone_home'])

            database_dict['add'].append(entry_dict)

        for tag in self.entries_to_remove.keys():
            entry = self.entries_to_remove[tag]
            database_dict.setdefault('remove', [])

            entry_dict = dataclasses.asdict(entry)

            entry_dict['value'] = entry_dict['tag']
            del entry_dict['tag']

            entry_dict['name'] = entry_dict['description']
            del entry_dict['description']

            if 'default' in entry_dict:
                entry_dict['default'] = str(entry_dict['default'])
            if 'phone_home' in entry_dict:
                entry_dict['phone_home'] = str(entry_dict['phone_home'])

            database_dict['remove'].append(entry_dict)

        return update_dict





# {
#   "data_tables": {
#     "testbeta": {
#       "add": [
#         {
#           "value": "newvalue",
#           "path": "newvalue.txt"
#         },
#         {
#           "value": "newvalue2",
#           "path": "newvalue2.txt"
#         }
#       ],
#       "remove": [
#         {
#           "value": "newvalue",
#           "path": "newvalue.txt"
#         }
#       ]
#     }
#   }
# }
        


# --context configure
# --mode add
# --mode update
# --mode remove

# --context manage
# --mode download
# --mode delete

def manage_database(args, params, locfile):
    entry = locfile.get_entry(args.tag)
    if entry is None:
        sys.exit(f'invalid database tag: {args.tag}')

        


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--context', required=True)
    parser.add_argument('--mode', required=True)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--output_file', required=True)
    parser.add_argument('--description')
    parser.add_argument('--default', default=0)
    parser.add_argument('--manifest')
    parser.add_argument('--local_cache_dir')
    parser.add_argument('--tool_cache_dir')
    parser.add_argument('--phone_home', default=0)

    args = parser.parse_args()

    with open(args.output_file) as f:
        params = json.load(f)

    with open('/var/tmp/dump.json', 'w') as ofh:
        print(json.dumps(params, indent=4), file=ofh)

    locfile_dir = params['param_dict']['__tool_directory__']
    locfile_path = os.path.join(locfile_dir, '..', 'test-data', 'ncbi_fcs_gx_databases.loc')
    locfile = LocFile(locfile_path)

    if args.context == 'configure':
        locfile.configure_database(args, params)
    elif args.context == 'manage':
        manage_database(args, params, locfile)
        



################################################################################

#    loc_file_dir = params[param_dict]['GALAXY_DATA_INDEX_DIR']

#    loc_file = LocFile(loc_file_path)


# ## NCBI FCS GX Databases
# # 
# #tag	description	default	manifest	local_cache_dir	tool_cache_dir
# all	Complete GX database	1	https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest	/scratch/rico/data	/dev/shm/gxdb
# test-only	Testing GX database	0	https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest	/scratch/rico/data	/dev/shm/gxdb








#    with open('/var/tmp/rico.json', 'w') as ofh:
#        print(json.dumps(params, indent=4), file=ofh)


if __name__ == '__main__':
    main()


# {
#   "param_dict": {
#     "chromInfo": "/scratch/rico/galaxy/tool-data/shared/ucsc/chrom/?.len",
#     "__datatypes_config__": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/registry.xml",
#     "GALAXY_DATATYPES_CONF_FILE": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/registry.xml",
#     "__history_id__": "359fec87b46b6eee",
#     "__galaxy_url__": "http://localhost:8080",
#     "database": "all",
#     "dbkey": "?",
#     "__input_ext": "input",
#     "__user__": "galaxy.model:SafeStringWrapper__galaxy.model.User__NoneType__NotImplementedType__Number__SafeStringWrapper__ToolParameterValueWrapper__bool__bytearray__ellipsis",
#     "__user_id__": "1",
#     "userId": "1",
#     "__user_email__": "planemo@galaxyproject.org",
#     "userEmail": "planemo@galaxyproject.org",
#     "__user_name__": "planemo",
#     "output_file": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/outputs/dataset_1a78f4d2-de5f-46fb-af98-8a7094de72cf.dat",
#     "__tool_directory__": "/scratch/rico/tool/data_managers/data_manager_ncbi_fcs_gx_database_downloader/data_manager",
#     "__get_data_table_entry__": "<function ToolEvaluator.__populate_non_job_params.<locals>.get_data_table_entry at 0x7f42897bc1e0>",
#     "__local_working_directory__": "/tmp/tmp3wuq6ckm/job_working_directory/000/1",
#     "__app__": "galaxy.app:UniverseApplication",
#     "__new_file_path__": "/tmp/tmp3wuq6ckm/tmp",
#     "__tool_data_path__": "/scratch/rico/galaxy/tool-data",
#     "GALAXY_DATA_INDEX_DIR": "/scratch/rico/galaxy/tool-data",
#     "__root_dir__": "/scratch/rico/galaxy",
#     "GALAXY_ROOT_DIR": "/scratch/rico/galaxy",
#     "__admin_users__": "planemo@galaxyproject.org,test@bx.psu.edu",
#     "input": "<function ToolEvaluator.build_param_dict.<locals>.input at 0x7f428994c0d0>"
#   },
#   "output_data": [
#     {
#       "out_data_name": "output_file",
#       "ext": "data_manager_json",
#       "dataset_id": 1,
#       "hda_id": 1,
#       "file_name": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/outputs/dataset_1a78f4d2-de5f-46fb-af98-8a7094de72cf.dat",
#       "extra_files_path": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/outputs/dataset_1a78f4d2-de5f-46fb-af98-8a7094de72cf_files"
#     }
#   ],
#   "job_config": {
#     "GALAXY_DATATYPES_CONF_FILE": "/tmp/tmp3wuq6ckm/job_working_directory/000/1/registry.xml",
#     "GALAXY_ROOT_DIR": "/scratch/rico/galaxy",
#     "TOOL_PROVIDED_JOB_METADATA_FILE": "galaxy.json"
#   }
# }
