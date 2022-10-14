import argparse
import errno
import hashlib
import json
import os
import subprocess as sp
import sys
import tarfile
from datetime import datetime
from pathlib import Path
from alive_progress import alive_bar
import requests
import bakta.constants as bc
import bakta.utils as bu


class GetBaktaDatabaseInfo:
    """
    Extract bakta database information to make a json file for data_manager
    """

    def __init__(self,
                 data_table_name="bakta_database",
                 output_path=Path.cwd().joinpath("db"),
                 db_version="latest"):
        self.bakta_table_list = None
        self.db_url = None
        self.data_table_entry = None
        self.data_table_name = data_table_name
        self.output_path = output_path
        self.db_version = db_version

    def get_data_table_format(self):
        """
        Build a data table format for galaxy
        using the bakta database information
        @str database_value: string of the database name
        @str database_date: string of the database date of build (YY-M-D)
        @str database_bakta_version: string of the version of bakta tool
         to apply a filter on version compatibility
        @str database_path: string of the database path
        for the database location
        return: a data table formatted for json output
        """
        self.data_table_entry = {
            "data_tables": {
                self.data_table_name: []
            }
        }
        return self.data_table_entry

    def fetch_db_versions(self, db_version="latest"):
        """
        Use method from bakta tool to extract database info
        db_version: a string of the version number
        in the galaxy wrapper list or just latest
        return: info for the select or the latest bakta db version
        """
        try:
            with requests.get(bc.DB_VERSIONS_URL) as resp:
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
                max(db_date_list)
                filtered_version = next(item for item in versions
                                        if max(db_date_list))
            elif db_version == "test":
                filtered_version = {"date": "date_test",
                                    "major": 0,
                                    "minor": 0,
                                    "doi": "10.5281/zenodo.7197299",
                                    "record": "7197299",
                                    "md5": "8b0250c17078742fc12207d4efb0fc1a",
                                    "software-min": {"major": 0, "minor": 0}}

            else:
                major_version = int(db_version.split(sep=".")[0])
                minor_version = int(db_version.split(sep=".")[1])
                try:
                    filtered_version = next(
                        item for item in versions
                        if item["major"] == major_version
                        and item["minor"] == minor_version)
                except StopIteration:
                    print("No available version detected in the list")
                    filtered_version = None
            self.db_url = f"https://zenodo.org/record/" \
                          f"{filtered_version['record']}/files/db.tar.gz"
            self.db_version = db_version
            return filtered_version

    def get_data_manager(self, bakta_database_info):
        self.bakta_table_list = self.get_data_table_format()
        data_info = dict(value=f"bakta_{bakta_database_info['major']}."
                               f"{bakta_database_info['minor']}",
                         dbkey=bakta_database_info['date'],
                         database_record=bakta_database_info['record'],
                         bakta_version=float(
                             f"{bakta_database_info['software-min']['major']}."
                             f"{bakta_database_info['software-min']['minor']}"
                         ), path=self.output_path.as_posix())
        toto = self.bakta_table_list["data_tables"][self.data_table_name]
        toto.append(data_info)
        return self.bakta_table_list


class InstallBaktaDatabase(GetBaktaDatabaseInfo):
    """
    Download the bakta database,
    check md5 sum,
    untar the download db and update for the amrfinderplus database
    """

    def __init__(self,
                 output_path=Path.cwd().joinpath("db"),
                 tarball_path=Path.cwd().joinpath("db.tar.gz")
                 ):
        super().__init__()
        self.md5 = None
        self.output_path = output_path
        self.tarball_path = tarball_path
        bu.test_dependency(bu.DEPENDENCY_AMRFINDERPLUS)

    def download(self):
        try:
            with self.tarball_path.open('wb') as fh_out, \
                    requests.get(self.db_url, stream=True) as resp:
                total_length = resp.headers.get('content-length')
                if total_length is None:  # no content length header
                    with alive_bar() as bar:
                        for data in resp.iter_content(chunk_size=1024 * 1024):
                            fh_out.write(data)
                            bar()
                else:
                    total_length = int(int(total_length) / 1024)
                    with alive_bar(total=total_length) as bar:
                        for data in resp.iter_content(chunk_size=1024 * 1024):
                            fh_out.write(data)
                            bar(incr=len(data) / 1024)
            print(f'Download bakta database {self.db_version}')
        except IOError:
            print(f'ERROR: Could not download file from Zenodo!'
                  f' url={self.db_url}, path={self.output_path}')

    def untar(self):
        try:
            with self.tarball_path.open('rb') as fh_in, \
                    tarfile.open(fileobj=fh_in, mode='r:gz') as tar_file:
                tar_file.extractall(path=str(self.output_path))
                print('Untar the database')
        except OSError:
            sys.exit(f'ERROR: Could not extract {self.tarball_path} '
                     f'to {self.output_path}')

    def calc_md5_sum(self, buffer_size=1048576):
        self.md5 = self.fetch_db_versions(db_version=self.db_version)["md5"]
        md5 = hashlib.md5()
        with self.tarball_path.open('rb') as fh:
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

    def update_amrfinderplus_db(self):
        amrfinderplus_db_path = f"{self.output_path}/db/amrfinderplus-db"
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
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            print(f"ERROR: AMRFinderPlus failed! "
                  f"command: 'amrfinder_update --force_update"
                  f" --database {amrfinderplus_db_path}'")
        else:
            print("AMRFinderPlus database download")


def parse_arguments():
    # parse options and arguments
    arg_parser = argparse.ArgumentParser("data_manager_json")
    arg_parser.add_argument("data_manager_json")
    arg_parser.add_argument("-d", "--database_version",
                            help='Select the database version '
                                 '(major and minor eg. 4.0),'
                                 'default is the latest version',
                            default="latest",
                            required=True)
    arg_parser.add_argument("-o",
                            "--output",
                            help='output path '
                                 '(default = current working directory)',
                            default=Path.cwd(), required=False)
    return arg_parser.parse_args()


def main():
    all_args = parse_arguments()
    with open(all_args.data_manager_json) as fh:
        data_manager_input = json.load(fh)
        target_dir = data_manager_input['output_data'][0]['extra_files_path']
        try:
            os.mkdir(target_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(target_dir):
                pass
            else:
                raise
    bakta_upload = InstallBaktaDatabase()
    bakta_db = bakta_upload.fetch_db_versions(
        db_version=all_args.database_version)
    bakta_upload.output_path = all_args.output
    bakta_data_manager = bakta_upload.get_data_manager(bakta_db)
    bakta_upload.download()
    bakta_upload.calc_md5_sum()
    bakta_upload.untar()
    bakta_upload.update_amrfinderplus_db()
    with open(all_args.data_manager_json, 'w') as fh:
        json.dump(bakta_data_manager, fh, sort_keys=True)


if __name__ == '__main__':
    main()
