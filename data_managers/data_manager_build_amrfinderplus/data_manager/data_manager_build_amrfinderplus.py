import argparse
import json
import os
import re
import subprocess as sp
from datetime import datetime
from pathlib import Path


class GetDataManager:

    def __init__(self):
        self.data_table_name = "amrfinderplus_database"
        self._db_name = "amrfinderplus-db"
        self._db_path = Path().absolute()
        self._today = datetime.now().strftime("%Y-%m-%d_%H:%M")

    def get_data_table_format(self):
        """
        Skeleton of a data_table format
        return: a data table formatted for json output
        """
        self.data_table_entry = {
            "data_tables": {
                self.data_table_name: {}
            }
        }
        return self.data_table_entry

    def get_data_manager(self, amrfinderplus_version):
        self.amrfinderplus_table_list = self.get_data_table_format()

        data_info = dict(value=self._today,
                         name=amrfinderplus_version,
                         path=self._db_name)
        self.amrfinderplus_table_list["data_tables"][self.data_table_name] = [data_info]
        return self.amrfinderplus_table_list

    def update_amrfinderplus_db(self):
        amrfinderplus_db_path = Path(self._db_path).joinpath(self._db_name)
        cmd = [
            'amrfinder_update',
            '--database', str(amrfinderplus_db_path),
            '--force_update'
        ]
        print(cmd)
        proc = sp.run(
            cmd,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            print(
                f"ERROR: AMRFinderPlus failed! command: 'amrfinder_update --force_update --database {amrfinderplus_db_path}', error code: {proc.returncode}")
        else:
            return amrfinderplus_db_path

    def get_amrfinderplus_version(self, amrfinderplus_path):
        version_file = Path(f'{amrfinderplus_path}/latest/version.txt')
        with open(version_file, "r") as version:
            version_value = version.read()
        version_value = re.sub("\n", "", version_value)
        return version_value

    def parse_arguments(self):
        # parse options and arguments
        arg_parser = argparse.ArgumentParser()
        arg_parser.add_argument("data_manager_json")
        arg_parser.add_argument("-t", "--test", action='store_true',
                                help="option to test the script with an lighted database")
        return arg_parser.parse_args()

    def read_json_input_file(self, json_file_path):
        with open(json_file_path) as fh:
            params = json.load(fh)
        target_dir = params['output_data'][0]['extra_files_path']
        os.makedirs(target_dir)
        return Path(target_dir)

    def write_json_infos(self, json_file_path, data_manager_infos):
        with open(json_file_path, 'w') as fh:
            json.dump(data_manager_infos, fh, sort_keys=True)


def main():
    # init the class
    amrfinderplus_download = GetDataManager()
    # import the arguments
    all_args = amrfinderplus_download.parse_arguments()
    # read the json input from galaxy to define the db path
    path_to_download = amrfinderplus_download.read_json_input_file(json_file_path=all_args.data_manager_json)
    # change the path to th json information
    amrfinderplus_download._db_path = path_to_download
    # download the last amrfinderplus database
    amrfinder_output = amrfinderplus_download.update_amrfinderplus_db()
    # extract the version number of the database
    amrfinder_version = amrfinderplus_download.get_amrfinderplus_version(amrfinder_output)
    # make a dic with database information
    amrfinderplus_json_output = amrfinderplus_download.get_data_manager(amrfinder_version)
    amrfinderplus_download.write_json_infos(json_file_path=all_args.data_manager_json,
                                            data_manager_infos=amrfinderplus_json_output)


if __name__ == '__main__':
    main()
