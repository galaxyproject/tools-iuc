import argparse
import json
import os
import tarfile
from datetime import datetime
from pathlib import Path

import requests


class GetPlasmidfinderDataManager:
    """
    Create the json file with database information for galaxy data manager
    """

    def __init__(self,
                 plasmidfinder_database="plasmidfinder_database",
                 db_name="plasmidfinder-db",
                 plasmidfinder_version="latest"):
        self.data_table_name = plasmidfinder_database
        self._db_name = db_name
        self._plasmidfinder_version = plasmidfinder_version
        self._plasmidfinder_date_version = None
        self.data_table_entry = None
        self.plasmidfinder_table_list = None
        self._commit_number = None

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

    def get_data_manager(self):
        """
        Create the empty data table format and add all the information into
        Commit number is added if latest is required instead of version number
        return: The data table with database information
        """
        self.plasmidfinder_table_list = self.get_data_table_format()
        if self._plasmidfinder_version == "latest":
            version_value = self._commit_number
        else:
            version_value = self._plasmidfinder_version
        plasmidfinder_value = f"plasmidfinder_{self._commit_number}" \
                              f"_{self._plasmidfinder_date_version}"
        plasmidfinder_name = f"{version_value}" \
                             f"_{self._plasmidfinder_date_version}"
        data_info = dict(value=plasmidfinder_value,
                         name=plasmidfinder_name,
                         date=self._plasmidfinder_date_version,
                         path=self._db_name)
        self.plasmidfinder_table_list["data_tables"][self.data_table_name] = [data_info]
        return self.plasmidfinder_table_list


class DownloadPlasmidfinderDatabase(GetPlasmidfinderDataManager):
    """
    Download the plasmidfinder database from the bitbucket repository.
    Build the data manager info for galaxy
    """

    def __init__(self,
                 output_dir=Path.cwd(),
                 plasmidfinder_url="https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/",
                 db_name="plasmidfinder-db",
                 db_tmp="tmp_database",
                 plasmidfinder_version="latest",
                 json_file_path=None,
                 date_version=datetime.now().strftime("%Y-%m-%d")):

        super().__init__()
        self.json_file_path = json_file_path
        self._output_dir = output_dir
        self._plasmidfinder_url = plasmidfinder_url
        self._temporary_folder = db_tmp
        self._db_name = db_name
        self._db_name_tar = f'{db_name}.gz'
        self._plasmidfinder_version = plasmidfinder_version
        self._plasmidfinder_date_version = date_version
        self._commit_number = None

    def extract_db_commit(self, request_header, title_name="content-disposition"):
        """
        Extract the commit if to add the information as identifier of the download
        @request_header: a request object obtained from requests.get()
        @title_name: the tag to search in the header of the requests object
        return: the value of the commit
        """
        db_info = request_header.headers[title_name]
        commit_number = db_info.split("-")[2].split(".")[0]
        return commit_number

    def untar_files(self, file_path: Path, extracted_path_output: Path):
        """
        untar the download archive
        @file_path: input path of the tar.gz file
        @extracted_path_output: output path of the extract folder
        return: the path of the output
        """
        try:
            with file_path.open('rb') as fh_in, \
                    tarfile.open(fileobj=fh_in, mode='r:gz') as tar_file:
                tar_file.extractall(path=extracted_path_output)
                print(f'Untar the database in {extracted_path_output}')
                return extracted_path_output
        except OSError:
            os.sys.exit(f'ERROR: Could not extract {file_path}')

    def choose_db_version(self):
        """
        Update the url link depending of the version choosen by user.
        This method could be upgraded simply by adding the new versions
        """
        if self._plasmidfinder_version == "latest":
            self._plasmidfinder_url = f"{self._plasmidfinder_url}master.gz"
        elif self._plasmidfinder_version == "2.1":
            self._plasmidfinder_url = f"{self._plasmidfinder_url}1307168.gz"

    def download_database(self):
        """
        Download the plasmidfinder database using requests lib
        Make the directory and temporary directory for download
        Untar the download files
        """
        self._output_dir = Path(self._output_dir)
        self.choose_db_version()
        try:
            request_info = requests.get(self._plasmidfinder_url)
            request_info.raise_for_status()
            self._commit_number = self.extract_db_commit(request_info)
            output_tar_path = self._output_dir.joinpath(self._temporary_folder)
            output_tar_path_file = output_tar_path.joinpath(self._db_name_tar)
            output_path = self._output_dir.joinpath(self._db_name)
            os.makedirs(output_tar_path)
            os.makedirs(output_path)
            with open(output_tar_path_file, 'wb') as output_dir:
                output_dir.write(request_info.content)
            untar_output = self.untar_files(file_path=output_tar_path_file, extracted_path_output=output_tar_path.joinpath(self._db_name))

            self.moove_download_files(older_path=untar_output, new_path=output_path)
        except requests.exceptions.HTTPError as http_error:
            print(f"Requests Error: {http_error}")
            print(f"Fail to import Plasmidfinder database from {self._plasmidfinder_url}")

    def moove_download_files(self, older_path, new_path, expression_search="*fsa"):
        """
        Clean downloaded data by mooving fasta files in the final folder
        @older_path: previous path where the files are located
        @new_path: final path where files will be mooved
        @expression_search: keep only file with this expression
        """
        fasta_files = Path(older_path).rglob(expression_search)
        file_list_paths = [file for file in fasta_files if file.is_file()]
        [self.keep_filename(pathname=path, output_path=new_path) for path in file_list_paths]

    def keep_filename(self, pathname, output_path):
        """
        Moove files
        @pathname: previous path
        @output_path: final path
        """
        Path.replace(pathname, output_path.joinpath(pathname.name))

    def read_json_input_file(self):
        """
        Import the json file
        """
        with open(self.json_file_path) as fh:
            params = json.load(fh)
        target_dir = params['output_data'][0]['extra_files_path']
        os.makedirs(target_dir)
        self._output_dir = target_dir

    def write_json_infos(self):
        """
        Write in the imported json file
        """
        with open(self.json_file_path, 'w') as fh:
            json.dump(self.get_data_manager(), fh, sort_keys=True)


def parse_arguments():
    """
    List of arguments provided by the user
    return: parsed arguments
    """
    # parse options and arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("data_manager_json",
                            help="json file from galaxy")
    arg_parser.add_argument("-v", "--db_version",
                            help="version of the plasmidfinder (latest or 2.1)")
    return arg_parser.parse_args()


def main():
    all_args = parse_arguments()
    plasmidfinder_download = DownloadPlasmidfinderDatabase(json_file_path=all_args.data_manager_json, plasmidfinder_version=all_args.db_version)
    plasmidfinder_download.read_json_input_file()
    plasmidfinder_download.download_database()
    plasmidfinder_download.write_json_infos()


if __name__ == '__main__':
    main()
