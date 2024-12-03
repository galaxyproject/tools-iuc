#!/usr/bin/env python

import argparse
import json
import os
import shutil
import tarfile
from datetime import datetime

import wget

DB_paths = {
    "mgnify_v5_lsu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_lsu-20200130.tar.gz",
    "mgnify_v5_ssu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_ssu-20200130.tar.gz",
    "mgnify_v5_its_unite": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/UNITE-20200214.tar.gz",
    "mgnify_v5_its_itsonedb": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/ITSoneDB-20200214.tar.gz",
    "mgnify_v6_lsu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-lsu/silva-lsu_138.1.tar.gz",
    "mgnify_v6_ssu": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/silva-ssu/silva-ssu_138.1.tar.gz",
    "mgnify_v6_its_unite": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/unite/unite_9.0.tar.gz",
    "mgnify_v6_its_itsonedb": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/itsonedb/itsonedb_1.141.tar.gz",
    "mgnify_v6_pr2": "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pr2/pr2_5.0.0.tar.gz",
    "test_lsu": "https://zenodo.org/record/8205348/files/test_lsu.tar.gz",
}

DB_names = {
    "mgnify_v5_lsu": "MGnify LSU (v5.0.7) - silva_lsu-20200130",
    "mgnify_v5_ssu": "MGnify SSU (v5.0.7) - silva_ssu-20200130",
    "mgnify_v5_its_unite": "MGnify ITS UNITE (v5.0.7) - UNITE-20200214",
    "mgnify_v5_its_itsonedb": "MGnify ITS ITSonedb (v5.0.7) - ITSoneDB-20200214",
    "mgnify_v6_lsu": "MGnify LSU (v6.0) - silva_lsu-20240702",
    "mgnify_v6_ssu": "MGnify SSU (v6.0) - silva_ssu-20240701",
    "mgnify_v6_its_unite": "MGnify ITS UNITE (v6.0) - UNITE-20240702",
    "mgnify_v6_its_itsonedb": "MGnify ITS ITSonedb (v6.0) - ITSoneDB-20240702",
    "mgnify_v6_pr2": "MGnify PR2 (v6.0) - PR2-20240702",
    "test_lsu": "Trimmed LSU Test DB",
}


def download_untar_store(url, tmp_path, dest_path):
    """
    Download a tar.gz file containing one folder,
    extract that folder and move the content inside dest_path
    """

    extract_path = os.path.join(tmp_path, "extract")

    os.makedirs(tmp_path, exist_ok=True)

    # download data
    filename = wget.download(url, out=tmp_path)
    tarfile_path = os.path.join(tmp_path, filename)
    tar = tarfile.open(tarfile_path)
    tar.extractall(extract_path)

    print(f"Content of folder: {extract_path}", os.listdir(extract_path))

    # case for mapseq v6: all DB files are directly in the tar.gz file
    # remove the VERSION.txt file since the tool can only handle on .txt file in the DB
    if len(list(os.listdir(extract_path))) > 1:
        print(f"Found multiple files in {extract_path}. Copy the content.")
        print(f"Copy data to {dest_path}")
        version_file_path = os.path.join(extract_path, "VERSION.txt")
        os.remove(version_file_path)
        shutil.copytree(extract_path, dest_path)
        print("Done !")

    # case for mapseq v5: all files are in a subfolder in the tar.gz file
    else:
        print(f"Found a folder in {extract_path}. Copy the content of the folder.")
        for folder in os.listdir(extract_path):
            folder_path = os.path.join(extract_path, folder)

            print(f"Copy data to {dest_path}")
            shutil.copytree(folder_path, dest_path)
            print("Done !")

    shutil.rmtree(tmp_path)


def main():
    # Parse Command Line
    parser = argparse.ArgumentParser(description="Create data manager JSON.")
    parser.add_argument("--out", dest="output", action="store", help="JSON filename")
    parser.add_argument(
        "--version", dest="version", action="store", help="Version of the DB"
    )
    parser.add_argument(
        "--database-type", dest="db_type", action="store", help="Db type"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="option to test the script with an lighted database",
    )

    args = parser.parse_args()

    # the output file of a DM is a json containing args that can be used by the DM
    # most tools mainly use these args to find the extra_files_path for the DM, which can be used
    # to store the DB data
    with open(args.output) as fh:
        params = json.load(fh)

    workdir = params["output_data"][0]["extra_files_path"]
    os.mkdir(workdir)

    time = datetime.utcnow().strftime("%Y-%m-%d")
    db_value = f"{args.db_type}_from_{time}"

    # output paths
    db_path = os.path.join(workdir, db_value)
    tmp_path = os.path.join(workdir, "tmp")

    # create DB
    if args.test:
        url = DB_paths["test_lsu"]
    else:
        url = DB_paths[args.db_type]

    # download data
    download_untar_store(url, tmp_path, db_path)

    db_name = DB_names[args.db_type]
    # Update Data Manager JSON and write to file
    data_manager_entry = {
        "data_tables": {
            "mapseq_db": {
                "value": db_value,
                "name": f"{db_name} downloaded at {time}",
                "version": args.version,
                "path": db_path,
            }
        }
    }

    with open(os.path.join(args.output), "w+") as fh:
        json.dump(data_manager_entry, fh, sort_keys=True)


if __name__ == "__main__":
    main()
