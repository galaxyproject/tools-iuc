import argparse
import hashlib
import json
import os
import sys
import tarfile
from datetime import datetime
from pathlib import Path


import requests


class GetBaktaDatabaseInfo:
    """
    Extract bakta database information to make a json file for data_manager
    """

    def __init__(self,
                 data_table_name="bakta_database",
                 db_name=Path.cwd().joinpath("db"),
                 db_version="latest",
                 test_mode=False):
        self.bakta_table_list = None
        self.db_url = None
        self.data_table_entry = None
        self.data_table_name = data_table_name
        self.db_name = db_name
        self.db_version = db_version
        self.DB_VERSIONS_URL = 'https://raw.githubusercontent.com/oschwengers/bakta/master/db-versions.json'
        self.DB_TEST_URL = 'https://zenodo.org/record/7360542/files/db-versions.json'
        self.test_mode = test_mode

    def get_data_table_format(self):
        """
        Skeleton of a data_table format
        return: a data table formated for json output
        """
        self.data_table_entry = {
            "data_tables": {
                self.data_table_name: {}
            }
        }
        return self.data_table_entry

    def fetch_db_versions(self, db_version="latest"):
        """
        List bakta database info related to the db_version selected
        """
        if self.test_mode is True:
            self.DB_VERSIONS_URL = self.DB_TEST_URL
        try:
            with requests.get(self.DB_VERSIONS_URL) as resp:
                versions = json.loads(resp.content)
        except IOError as e:
            print(e, file=sys.stderr)
            raise e
        else:
            if db_version == "latest":
                db_date_list = []
                for db_dic in versions:
                    db_date_list.append(datetime.strptime(db_dic["date"],
                                                          '%Y-%m-%d').date())
                filtered_version = max(versions, key=lambda x: x['date'])
            else:
                filtered_version = None
                for item in versions:
                    if '{0}.{1}'.format(item["major"], item["minor"]) == db_version:
                        filtered_version = item
                        break
                if filtered_version is None:
                    print("No matching version detected in the list")
            if filtered_version is not None:
                self.db_url = f"https://zenodo.org/record/" \
                              f"{filtered_version['record']}/files/db.tar.gz"
                self.db_version = db_version
                return filtered_version

    def get_data_manager(self, bakta_database_info):
        self.bakta_table_list = self.get_data_table_format()
        bakta_name = f"V{bakta_database_info['major']}." \
                     f"{bakta_database_info['minor']}_" \
                     f"{bakta_database_info['date']}"
        tool_version = str(f"{bakta_database_info['software-min']['major']}."
                           f"{bakta_database_info['software-min']['minor']}")
        data_info = dict(value=bakta_name,
                         dbkey=bakta_database_info['record'],
                         bakta_version=tool_version,
                         path="db")
        self.bakta_table_list["data_tables"][self.data_table_name] = [data_info]
        return self.bakta_table_list


class InstallBaktaDatabase(GetBaktaDatabaseInfo):
    """
    Download the bakta database,
    check md5 sum,
    untar the download db and update for the amrfinderplus database
    """

    def __init__(self,
                 db_dir=Path.cwd(),
                 db_name="bakta",
                 tarball_name="db.tar.gz",
                 test_mode=False):
        super().__init__()
        self.md5 = None
        self.db_dir = db_dir
        self.db_name = db_name
        self.tarball_name = tarball_name
        self.tarball_path = None
        self.test_mode = test_mode

    def download(self):
        self.db_name = f'{self.db_name}_{self.db_version}'
        bakta_path = Path(self.db_dir).joinpath(self.tarball_name)
        try:
            with bakta_path.open('wb') as fh_out, \
                    requests.get(self.db_url, stream=True) as resp:
                total_length = resp.headers.get('content-length')
                if total_length is None:  # no content length header
                    for data in resp.iter_content(chunk_size=1024 * 1024):
                        fh_out.write(data)
                else:
                    for data in resp.iter_content(chunk_size=1024 * 1024):
                        fh_out.write(data)
            print(f'Download bakta database {self.db_version}')
            self.tarball_path = bakta_path
        except IOError:
            print(f'ERROR: Could not download file from Zenodo!'
                  f' url={self.db_url}, path={self.tarball_name}')

    def untar(self):
        db_path = Path(self.db_dir).as_posix()
        try:
            with self.tarball_path.open('rb') as fh_in, \
                    tarfile.open(fileobj=fh_in, mode='r:gz') as tar_file:
                tar_file.extractall(path=db_path)
                print(f'Untar the database in {db_path}')
                return db_path
        except OSError:
            sys.exit(f'ERROR: Could not extract {self.tarball_name} '
                     f'to {self.db_name}')

    def calc_md5_sum(self, buffer_size=1048576):
        tarball_path = Path(self.db_dir).joinpath(self.tarball_name)
        self.md5 = self.fetch_db_versions(db_version=self.db_version)["md5"]
        md5 = hashlib.md5()
        with tarball_path.open('rb') as fh:
            data = fh.read(buffer_size)
            while data:
                md5.update(data)
                data = fh.read(buffer_size)
        if md5.hexdigest() == self.md5:
            print('\t...md5 control database OK')
        else:
            print(f"Error: corrupt database file! "
                  f"calculated md5 = {md5.hexdigest()}"
                  f" different from {self.md5} ")


"""
This is the method to download the amrfinderplus database need by bakta.
Deprecated to use the amrfinderplus data_manager
    def update_amrfinderplus_db(self):
        amrfinderplus_db_path = f"{self.db_dir}/{self.db_name}/db/amrfinderplus-db"
        if self.db_version == "test":
            cmd = [
                'amrfinder_update',
                '--database', str(amrfinderplus_db_path),
                '--force_update',
                '--help'
            ]
        else:
            cmd = [
                'amrfinder_update',
                '--database', str(amrfinderplus_db_path),
                '--force_update'
            ]
        proc = sp.run(
            cmd,
            universal_newlines=True
        )
        if proc.returncode != 0:
            print(f"ERROR: AMRFinderPlus failed! "
                  f"command: 'amrfinder_update --force_update"
                  f" --database {amrfinderplus_db_path}'")
        else:
            print("AMRFinderPlus database download")
"""


def parse_arguments():
    # parse options and arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("data_manager_json")
    arg_parser.add_argument("-d", "--database_version",
                            help='Select the database version '
                                 '(major and minor eg. 4.0),'
                                 'default is the latest version',
                            default="latest",
                            required=True)
    arg_parser.add_argument("-t", "--test", action='store_true',
                            help="option to test the script with an empty database")
    return arg_parser.parse_args()


def main():
    all_args = parse_arguments()
    with open(all_args.data_manager_json) as fh:
        params = json.load(fh)
    target_dir = params['output_data'][0]['extra_files_path']
    os.makedirs(target_dir)
    # init the class to download bakta db
    bakta_upload = InstallBaktaDatabase(test_mode=all_args.test)
    bakta_db = bakta_upload.fetch_db_versions(db_version=all_args.database_version)
    # update the path for galaxy
    bakta_upload.db_dir = target_dir
    # download the database
    bakta_upload.download()
    # check md5 sum
    bakta_upload.calc_md5_sum()
    # untar db
    bakta_upload.untar()
    # make the data_manager metadata
    bakta_data_manager = bakta_upload.get_data_manager(bakta_database_info=bakta_db)
    with open(all_args.data_manager_json, 'w') as fh:
        json.dump(bakta_data_manager, fh, sort_keys=True)


if __name__ == '__main__':
    main()
