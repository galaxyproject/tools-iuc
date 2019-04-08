#!/usr/bin/env python
# Data manager for reference data for the QIIME Galaxy tools

import argparse
import ftplib
import json
import os
import tarfile
import zipfile

import requests


protocol = {
    "unite": "http",
    "greengenes": "ftp",
    "silva": "http",
    "img": "ftp"
}
baseUrl = {
    "unite": "http://unite.ut.ee/sh_files/sh_qiime_release_",
    "greengenes": "greengenes.microbio.me",
    "silva": "http://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_",
    "img": "ftp.microbio.me"
}
ftp_dir = {
    "greengenes": "/greengenes_release/gg_",
    "img": ""
}
ftp_file_prefix = {
    "greengenes": "gg_",
    "img": ""
}
ftp_file_suffix = {
    "greengenes": "_otus",
    "img": ""
}
extension = {
    "unite": "zip",
    "greengenes": "tar.gz",
    "silva": {
        "104_release": "tgz",
        "108_release": "tgz",
        "108_release_curated": "tgz",
        "111_release": "tgz",
        "119_consensus_majority_taxonomy": "zip",
        "119_release": "zip",
        "119_release_aligned_rep_files": "tar.gz",
        "123_release": "zip",
        "128_release": "tgz"},
    "img": "tgz"
}
filetypes = ["rep_set", "rep_set_aligned", "taxonomy", "trees"]


# Utility functions for interacting with Galaxy JSON
def read_input_json(jsonfile):
    """Read the JSON supplied from the data manager tool

    Returns a tuple (param_dict,extra_files_path)

    'param_dict' is an arbitrary dictionary of parameters
    input into the tool; 'extra_files_path' is the path
    to a directory where output files must be put for the
    receiving data manager to pick them up.

    NB the directory pointed to by 'extra_files_path'
    doesn't exist initially, it is the job of the script
    to create it if necessary.

    """
    params = json.loads(open(jsonfile).read())
    return (params['param_dict'],
            params['output_data'][0]['extra_files_path'])


# Utility functions for creating data table dictionaries
#
# Example usage:
# >>> d = create_data_tables_dict()
# >>> add_data_table(d,'my_data')
# >>> add_data_table_entry(dict(dbkey='hg19',value='human'))
# >>> add_data_table_entry(dict(dbkey='mm9',value='mouse'))
# >>> print str(json.dumps(d))
def create_data_tables_dict():
    """Return a dictionary for storing data table information

    Returns a dictionary that can be used with 'add_data_table'
    and 'add_data_table_entry' to store information about a
    data table. It can be converted to JSON to be sent back to
    the data manager.

    """
    d = {}
    d['data_tables'] = {}
    return d


def add_data_table(d, table):
    """Add a data table to the data tables dictionary

    Creates a placeholder for a data table called 'table'.

    """
    d['data_tables'][table] = []


def add_data_table_entry(d, table, entry):
    """Add an entry to a data table

    Appends an entry to the data table 'table'. 'entry'
    should be a dictionary where the keys are the names of
    columns in the data table.

    Raises an exception if the named data table doesn't
    exist.

    """
    try:
        d['data_tables'][table].append(entry)
    except KeyError:
        raise Exception("add_data_table_entry: no table '%s'" % table)


def get_ftp_file(ftp, filename):
    """
    """
    try:
        ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)
    except Exception:
        print("Error")


def download_archive(db, version, ext):
    """

    """
    filepath = "%s_%s.%s" % (db, version, ext)
    if protocol[db] == "http":
        url = "%s%s.%s" % (baseUrl[db], version, ext)
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(filepath, "wb") as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)
    elif protocol[db] == "ftp":
        ftp = ftplib.FTP(baseUrl[db])
        ftp.login("anonymous", "ftplib-example-1")
        if db == "greengenes" and version == "13_8":
            ftp.cwd("%s%s" % (ftp_dir[db], "13_5"))
        else:
            ftp.cwd("%s%s" % (ftp_dir[db], version))
        filepath = "%s%s%s.%s" % (
            ftp_file_prefix[db],
            version,
            ftp_file_suffix[db],
            ext)
        get_ftp_file(ftp, filepath)
        ftp.quit()
    return filepath


def find_archive_content_path(archive_content_path):
    """
    """
    content = os.listdir(archive_content_path)
    archive_content = []
    for x in content:
        if not x.startswith(".") and not x.startswith("_"):
            archive_content.append(x)
    if len(archive_content) == 1:
        archive_content_path = os.path.join(
            archive_content_path,
            archive_content[0])
    return archive_content_path


def extract_archive(filepath, ext, db):
    """
    """
    archive_content_path = "tmp"
    if ext == "tar.gz" or ext == "tgz":
        tar = tarfile.open(filepath)
        tar.extractall(path=archive_content_path)
        tar.close()
        archive_content_path = find_archive_content_path(archive_content_path)
    elif ext == "zip":
        zip_ref = zipfile.ZipFile(filepath, 'r')
        zip_ref.extractall(archive_content_path)
        zip_ref.close()
        archive_content_path = find_archive_content_path(archive_content_path)
    return archive_content_path


def move_unite_files(archive_content_path, filename_prefix, name_prefix, data_tables, target_dir):
    """

    """
    archive_content = os.listdir(archive_content_path)
    for content in archive_content:
        content_filepath = os.path.join(archive_content_path, content)
        content_name_prefix = "%s - %s" % (name_prefix, content.split(".")[0])
        content_filename_prefix = "%s_%s" % (filename_prefix, content)
        if content.find("refs") != -1:
            move_file(
                content_filepath,
                content_filename_prefix,
                content_name_prefix,
                data_tables,
                os.path.join(target_dir, "rep_set"),
                "rep_set")
        elif content.find("taxonomy") != -1:
            move_file(
                content_filepath,
                content_filename_prefix,
                content_name_prefix,
                data_tables,
                os.path.join(target_dir, "taxonomy"),
                "taxonomy")


def move_file(input_filepath, filename, name, data_tables, target_dir, filetype):
    """
    """
    output_filepath = os.path.join(target_dir, filename)
    os.rename(input_filepath, output_filepath)
    add_data_table_entry(
        data_tables,
        "qiime_%s" % (filetype),
        dict(
            dbkey=filename,
            value=os.path.splitext(filename)[0],
            name=name,
            path=output_filepath))


def move_dir_content(input_path, filename_prefix, name_prefix, data_tables, target_dir, filetype):
    """
    """
    for content in os.listdir(input_path):
        if content.startswith("."):
            continue
        content_path = os.path.join(input_path, content)
        content_name_prefix = "%s - %s" % (name_prefix, content.split(".")[0])
        content_filename_prefix = "%s_%s" % (filename_prefix, content)
        if os.path.isdir(content_path):
            move_dir_content(
                content_path,
                content_filename_prefix,
                content_name_prefix,
                data_tables,
                target_dir,
                filetype)
        else:
            move_file(
                content_path,
                content_filename_prefix,
                content_name_prefix,
                data_tables,
                target_dir,
                filetype)


def move_files(archive_content_path, filename_prefix, name_prefix, data_tables, target_dir, db, version):
    """
    """
    for filetype in filetypes:
        if filetype == "rep_set_aligned":
            if db == "greengenes" and version == "12_10":
                continue
        filetype_target_dir = os.path.join(
            target_dir,
            filetype)
        filetype_path = os.path.join(
            archive_content_path,
            filetype)
        move_dir_content(
            filetype_path,
            filename_prefix,
            name_prefix,
            data_tables,
            filetype_target_dir,
            filetype)


def download_db(data_tables, db, version, target_dir):
    """Download QIIME database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      db: name of the database
      version: version of the database
      table_name: name of the table
      target_dir: directory to put copy or link to the data file

    """
    ext = extension[db]
    if db == "silva":
        ext = ext[version]

    print("Download archive")
    filepath = download_archive(db, version, ext)

    print("Extract archive %s" % filepath)
    archive_content_path = extract_archive(filepath, ext, db)

    print("Moving file from %s" % archive_content_path)
    filename_prefix = "%s_%s" % (db, version)
    name_prefix = "%s (%s)" % (db, version)
    if db == "greengenes" or db == "silva":
        move_files(
            archive_content_path,
            filename_prefix,
            name_prefix,
            data_tables,
            target_dir,
            db,
            version)
    elif db == "unite":
        move_unite_files(
            archive_content_path,
            filename_prefix,
            name_prefix,
            data_tables,
            target_dir)


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser(
        description='Download QIIME reference database')
    parser.add_argument('--database', help="Database name")
    parser.add_argument('--version', help="Database version")
    parser.add_argument('--jsonfile', help="Output JSON file")
    args = parser.parse_args()

    jsonfile = args.jsonfile

    # Read the input JSON
    params, target_dir = read_input_json(jsonfile)

    # Make the target directory
    print("Making %s" % target_dir)
    os.mkdir(target_dir, mode=0o755)
    os.mkdir(os.path.join(target_dir, "rep_set"), mode=0o755)
    os.mkdir(os.path.join(target_dir, "rep_set_aligned"), mode=0o755)
    os.mkdir(os.path.join(target_dir, "taxonomy"), mode=0o755)
    os.mkdir(os.path.join(target_dir, "trees"), mode=0o755)

    # Set up data tables dictionary
    data_tables = create_data_tables_dict()
    add_data_table(data_tables, "qiime_rep_set")
    add_data_table(data_tables, "qiime_rep_set_aligned")
    add_data_table(data_tables, "qiime_taxonomy")
    add_data_table(data_tables, "qiime_trees")

    # Fetch data from specified data sources
    download_db(
        data_tables,
        args.database,
        args.version,
        target_dir)

    # Write output JSON
    print("Outputting JSON")
    print(str(json.dumps(data_tables)))
    with open(jsonfile, 'w') as out:
        json.dump(data_tables, out)
    print("Done.")
