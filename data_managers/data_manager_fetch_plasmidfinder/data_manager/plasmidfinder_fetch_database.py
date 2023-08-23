import argparse
import json
import os
import time
from pathlib import Path


import git


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
                 plasmidfinder_url="https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master",
                 db_name="plasmidfinder-db",
                 db_tmp="tmp_database",
                 plasmidfinder_version="latest",
                 json_file_path=None,
                 date_version=None):

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

    def git_clone(self):
        git.Repo.clone_from(url=self._plasmidfinder_url, to_path=self._output_dir)
        self._plasmidfinder_repository = git.Repo(path=self._output_dir)

    def get_commit_number(self):
        sha = self._plasmidfinder_repository.head.commit.hexsha
        short_sha = self._plasmidfinder_repository.git.rev_parse(sha, short=7)
        self._commit_number = short_sha

    def get_commit_date(self):
        self._plasmidfinder_date_version = time.strftime("%Y_%m_%d", time.gmtime(self._plasmidfinder_repository.head.commit.committed_date))

    def download_database(self):
        """
        Download the plasmidfinder database using git lib
        Extract commit and commit date
        """
        self._output_dir = Path(self._output_dir)
        self.git_clone()
        if self._plasmidfinder_version != "latest":
            self._plasmidfinder_repository.git.checkout(self._plasmidfinder_version)
        self.get_commit_number()
        self.get_commit_date()

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