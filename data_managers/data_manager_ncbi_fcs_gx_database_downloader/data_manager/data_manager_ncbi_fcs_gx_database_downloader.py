#!/usr/bin/env python

import argparse
import json
import os
import sys
import subprocess
import dataclasses


# @dataclasses.dataclass
# class LocFileEntry:
#     tag: str
#     description: str
#     default: int = 0
#     manifest: str
#     local_cache_dir: str
#     tool_cache_dir: str

class LocFile:
    def __init__(self, pathname):
        self.pathname = pathname
        self.comments = []
        self.entries = {}
        self._load_file()

    def _load_file():
        with open(self.pathname) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    self.comments.append(line)
                else:
                    tag, description, default, manifest, local_cache_dir, tool_cache_dir = line.split('\t', 5)
                    default = int(default)
                    entry = LocFileEntry(tag, description, default, manifest, local_cache_dir, tool_cache_dir)
                    self.add_entry(entry)

    def add_entry(self, entry):
        self.entries[entry.tag] = entry
#        if entry.tag in self.entries:
#            # should this be an error??
#            self.update_entry(entry)
#        else:
#            self.entries[entry.tag] = entry

# python '$__tool_directory__/data_manager_ncbi_fcs_gx_database_downloader.py'
# --'${context.context_selector}'
# --'${context.mode.mode_selector}'
# #if $context.context_selector == 'configure'
#     #if $varExists('context.mode.tag')
#     --tag '${context.mode.tag}'
#     #else
#     --tag '${context.mode.database}'
#     #end if
# 
#     #if $context.mode_selector in ('add', 'update')
#     --description '${context.mode.description}'
#     --default '${context.mode.default}'
#     --manifest '${context.mode.manifest}'
#     --local_cache_dir '${context.mode.local_cache_dir}'
#     --tool_cache_dir '${context.mode.tool_cache_dir}'
#     #end if
# #else if $context.context_selector == 'manage'
#     --tag '${context.mode.database}'
# #end if
# --output_file '${output_file}'

# --configure --add
# --configure --update
# --configure --delete
# 
# --manage --download
# --manage --delete
        


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    parser.add_argument('--output_file', required=True)
    subparsers = parser.add_subparsers()

    parser_configure = subparsers.add_parser('--configure')
    parser_configure.add_argument('--description')
    parser_configure.add_argument('--default')
    parser_configure.add_argument('--manifest')
    parser_configure.add_argument('--local_cache_dir')
    parser_configure.add_argument('--tool_cache_dir')

    parser_manage = subparsers.add_parser('--manage')

    args = parser.parse_args()

    print (args)



#    with open(args.output_file) as f:
#    params = json.load(f)

################################################################################

#    loc_file_dir = params[param_dict]['GALAXY_DATA_INDEX_DIR']
#    loc_file_path = os.path.join(loc_file_dir, 'ncbi_fcs_gx_databases.loc')

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
