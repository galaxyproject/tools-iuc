import argparse
import json
import os
import subprocess as sp
from ftplib import FTP

from io import BytesIO
from pathlib import Path

import pandas as pd

class GetAmrFinderPlusDataManager:
    """
    Create the json file with database information for galaxy data manager
    """

    def __init__(self,
                 amrfinderplus_database="amrfinderplus_database",
                 db_name="amrfinderplus-db",
                 amrfinderplus_version="latest",
                 date_version=None):
        self.data_table_name = amrfinderplus_database
        self._db_name = db_name
        self._amrfinderplus_version = amrfinderplus_version
        self._amrfinderplus_date_version = date_version
        self.data_table_entry = None
        self.amrfinderplus_table_list = None

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
        return: The data table with database information
        """
        self.amrfinderplus_table_list = self.get_data_table_format()
        amrfinderplus_value = f"amrfinderplus_V{self._amrfinderplus_version}" \
                              f"_{self._amrfinderplus_date_version}"
        amrfinderplus_name = f"V{self._amrfinderplus_version}" \
                             f"-{self._amrfinderplus_date_version}"
        data_info = dict(value=amrfinderplus_value,
                         name=amrfinderplus_name,
                         path=self._db_name)
        self.amrfinderplus_table_list["data_tables"][self.data_table_name] = [data_info]
        return self.amrfinderplus_table_list


class DownloadAmrFinderPlusDatabase(GetAmrFinderPlusDataManager):
    """
    Download the amrfinderplus database from the ncbi.
    Make the database available with hmm and indexed files
    Build the data manager infos for galaxy
    """

    def __init__(self,
                 output_dir=Path.cwd(),
                 ncbi_url="ftp.ncbi.nlm.nih.gov",
                 ftp_login="anonymous",
                 ftp_password="anonymous",
                 amrfinderplus_database="amrfinderplus_database",
                 db_name="amrfinderplus-db",
                 amrfinderplus_version="latest",
                 json_file_path=None,
                 date_version=None,
                 amrfinderplus_db_path=None,
                 test_mode=False):

        super().__init__()
        self.json_file_path = json_file_path
        self._output_dir = output_dir
        self._ncbi_ftp_url = ncbi_url
        self._ncbi_database_path = "pathogen/Antimicrobial_resistance/AMRFinderPlus/database"
        self._login = ftp_login
        self._password = ftp_password
        self._amrfinderplus_database = amrfinderplus_database
        self._db_name = db_name
        self._amrfinderplus_version = amrfinderplus_version
        self._amrfinderplus_date_version = date_version
        self.species_list = None
        self.test_mode = test_mode
        self.amrfinderplus_db_path = amrfinderplus_db_path

    @staticmethod
    def subprocess_cmd(command, *args):
        """
        Method to call external tools with any parameters
        :param command: command name from the tool used (e.g. wget or makeblastdb)
        :param args: free number of argument need for the command tool (e.g. -r, -P ...)
        :return: launch the command line from the system
        """
        cmd = [command]
        [cmd.append(i) for i in args]
        proc = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
        if proc.returncode != 0:
            print(f'Error type {proc.returncode} with : \n {proc}')

    def download_amrfinderplus_db(self):
        """
        Download the amrfinderplus database from the ncbi ftp server
        """
        self.amrfinderplus_db_path = f'{self._output_dir}/{self._db_name}'
        os.makedirs(self.amrfinderplus_db_path)
        if self._amrfinderplus_version == 'latest':
            self.get_amrfinderplus_version()

        amrfinderplus_ftp_path = f"ftp://{self._login}:" \
                                 f"{self._password}@{self._ncbi_ftp_url}/" \
                                 f"{self._ncbi_database_path}/" \
                                 f"{self._amrfinderplus_version}/" \
                                 f"{self._amrfinderplus_date_version}"
        if self.test_mode is True:
            file_list = ["AMR_DNA-Escherichia", "version.txt", "taxgroup.tab", "database_format_version.txt"]
            output_option = "-O"
            for file in file_list:
                self.subprocess_cmd("wget",
                                    "-nd",
                                    "-np",
                                    "-r",
                                    f"{amrfinderplus_ftp_path}/{file}",
                                    output_option,
                                    f"{self.amrfinderplus_db_path}/{file}")
        else:
            output_option = "-P"
            self.subprocess_cmd("wget",
                                "-nd",
                                "-np",
                                "-r",
                                amrfinderplus_ftp_path,
                                output_option,
                                self.amrfinderplus_db_path)

    def make_hmm_profile(self):
        """
        Make the hmm profile using the AMR.LIB file previously download
        """
        hmm_file = Path(f"{self.amrfinderplus_db_path}/AMR.LIB")
        if Path.exists(hmm_file) and self.test_mode is False:
            self.subprocess_cmd("hmmpress", "-f", hmm_file)
        else:
            print("hmm_file file is missing to make hmm profiles")

    def extract_filelist_makeblast(self):
        """
        Extract le list of species which have file in the database
        return: a filtered species list of available species in the database
        """
        taxa_group_path = Path(f"{self.amrfinderplus_db_path}/taxgroup.tab")
        if Path.exists(taxa_group_path):
            taxa_table = pd.read_table(taxa_group_path)
            taxa_table.columns = ["taxgroup", "gpipe_taxgroup", "number_of_nucl_ref_genes"]
            taxa_df = taxa_table[taxa_table.number_of_nucl_ref_genes > 0].filter(items=["taxgroup"], axis=1)
            if self.test_mode is True:
                taxa_df = taxa_df[taxa_df.taxgroup == "Escherichia"].taxgroup
            else:
                taxa_df = taxa_df.taxgroup
            self.species_list = list(taxa_df)
        else:
            print("taxgroup.tab file is missing to list available species")

    def make_blastdb(self):
        """
        Index fasta file for blast
        """
        self.extract_filelist_makeblast()
        nucl_file_db_list = [f'{self.amrfinderplus_db_path}/AMR_DNA-{specie}' for specie in self.species_list]
        amr_dna = f'{self.amrfinderplus_db_path}/AMR_CDS'
        amr_prot = f'{self.amrfinderplus_db_path}/AMRProt'
        os.chdir(self.amrfinderplus_db_path)
        if Path(amr_dna).exists():
            nucl_file_db_list.append(amr_dna)
        else:
            print("No file AMR_CDS detected for indexing")
        if Path(amr_prot).exists():
            self.subprocess_cmd("makeblastdb", "-in", amr_prot, "-dbtype", "prot")
        else:
            print("No file AMRProt detected for indexing")
        [self.subprocess_cmd("makeblastdb", "-in", file, "-dbtype", "nucl") for file in nucl_file_db_list]

    def get_amrfinderplus_version(self, version_file="version.txt",
                                  database_version_file="database_format_version.txt"):
        """
        Check the version when latest if provided and update the number
        param version_file: name of the file containing version information
        param database_version_file: name of the file containing date version information
        """
        ftp = FTP(self._ncbi_ftp_url)
        ftp.login(self._login, self._password)
        ftp.cwd(f"{self._ncbi_database_path}/{self._amrfinderplus_version}")
        db_version = BytesIO()
        db_date_version = BytesIO()
        ftp.retrbinary(f'RETR {version_file}', db_version.write)
        ftp.retrbinary(f'RETR {database_version_file}', db_date_version.write)
        self._amrfinderplus_date_version = db_version.getvalue().decode("utf-8").splitlines()[0]
        self._amrfinderplus_version = '.'.join(
            db_date_version.getvalue().decode("utf-8").splitlines()[0].split(".")[:2])

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
    arg_parser.add_argument("--db_version", default="latest",
                            help="select the major version of the database (e.g. 3.10, 3.8), default is latest")
    arg_parser.add_argument("--db_date",
                            help="select the date into the database version (e.g. 2022-10-11.2)")
    arg_parser.add_argument("--test", action='store_true',
                            help="option to test the script with an lighted database")
    return arg_parser.parse_args()


def main():
    all_args = parse_arguments()
    amrfinderplus_download = DownloadAmrFinderPlusDatabase(amrfinderplus_version=all_args.db_version,
                                                           date_version=all_args.db_date,
                                                           json_file_path=all_args.data_manager_json,
                                                           test_mode=all_args.test)
    amrfinderplus_download.read_json_input_file()
    amrfinderplus_download.download_amrfinderplus_db()
    amrfinderplus_download.make_hmm_profile()
    amrfinderplus_download.make_blastdb()
    amrfinderplus_download.write_json_infos()


if __name__ == '__main__':
    main()
