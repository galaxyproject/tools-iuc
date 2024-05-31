#!/usr/bin/env python
#
# Data manager for reference data for the MetaPhlAn Galaxy tools
import argparse
import json
import subprocess
from datetime import date
from pathlib import Path


# Utility functions for interacting with Galaxy JSON
def read_input_json(json_fp):
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
    with open(json_fp) as fh:
        params = json.load(fh)
    return (params['param_dict'],
            Path(params['output_data'][0]['extra_files_path']))


# Utility functions for creating data table dictionaries
#
# Example usage:
# >>> d = create_data_tables_dict()
# >>> add_data_table(d,'my_data')
# >>> add_data_table_entry(dict(dbkey='hg19',value='human'))
# >>> add_data_table_entry(dict(dbkey='mm9',value='mouse'))
# >>> print(json.dumps(d))
def create_data_tables_dict():
    """Return a dictionary for storing data table information

    Returns a dictionary that can be used with 'add_data_table'
    and 'add_data_table_entry' to store information about a
    data table. It can be converted to JSON to be sent back to
    the data manager.

    """
    d = {
        'data_tables': {}
    }
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


def download_groot_db(data_tables, name, table_name, target_dp, identity, groot_version):
    """Download GROOT database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      name: name of the database to download
      table_name: name of the table
      target_dp: directory to put copy or link to the data file
      identity: identity threshold for GROOT
      groot_version: version of GROOT to use

    """
    # Define the target directory path
    db_dp = target_dp / Path(name)

    # Build the command string
    cmd = f"groot get -d %s -o %s --identity %s" % (name, db_dp, identity)

    # Execute the command
    subprocess.check_call(cmd, shell=True)

    # Add the data table entry
    add_data_table_entry(
        data_tables,
        table_name,
        dict(
            value='%s.%s-v%s' % (name, identity, groot_version),
            name='%s (%s percent identity)' % (name, identity),
            dbkey='%s-v%s' % (date.today().strftime("%d%m%Y"), groot_version),
            path=str(db_dp),
            db_version=groot_version
        )
    )


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser(description='Download and build Groot database')
    parser.add_argument('--database', help="Name of the database")
    parser.add_argument('--percentidentity', help="The identity threshold at which the database was clustered")
    parser.add_argument('--grootversion', help="Version of the Database")
    parser.add_argument('--json', help="Path to JSON file")
    args = parser.parse_args()
    print("args   : %s" % args)

    # Read the input JSON
    json_fp = Path(args.json)
    params, target_dp = read_input_json(json_fp)

    # Make the target directory
    print("Making %s" % target_dp)
    target_dp.mkdir(parents=True, exist_ok=True)

    # Set up data tables dictionary
    
    data_tables = create_data_tables_dict()
    add_data_table(data_tables, "groot_database_downloader")

    # Fetch data from specified data sources
    print("Download and build database")
    download_groot_db(
        data_tables,
        args.database,
        "groot_database_downloader",
        target_dp,
        args.percentidentity,
        args.grootversion)

    # Write output JSON
    print("Outputting JSON")
    with open(json_fp, 'w') as fh:
        json.dump(data_tables, fh, sort_keys=True)
    print("Done.")