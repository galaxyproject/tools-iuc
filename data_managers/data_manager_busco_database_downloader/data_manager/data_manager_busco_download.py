#!/usr/bin/env python
#
# Data manager for reference data for the 'BUSCO' Galaxy tools
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


def download_busco_db(data_tables, table_name, database, version, index, target_dp):
    """Download BUSCO database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      table_name: name of the table
      database: database to download (chocophlan or uniref)
      version: tool version
      target_dp: directory to put copy or link to the data file
    """
    db_target_dp = target_dp / Path(database)
    db_dp = db_target_dp / Path(database)
    # launch tool to get db
    cmd = "busco --download %s" % database
    subprocess.check_call(cmd, shell=True)
    # move db
    db_dp.rename(db_dp)
    # add details to data table
    add_data_table_entry(
        data_tables,
        table_name,
        dict(
            dbdate="%s" % date.today().strftime("%d%m%Y")  ,
            name="BUSCO database %s" % database,
            dbversion="db-version-%s-%s" % (version,index),
            path=str(db_dp)))


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser(description='Download BUSCO database')
    parser.add_argument('--database', help="Database name")
    parser.add_argument('--index', help="BUSCO database version")
    parser.add_argument('--version', help="BUSCO version")
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
    if args.database == "all":
        table_name = 'busco_all_database'
    elif args.database == "prokaryota":
        table_name = 'busco_prokaryota_database'
    elif args.database == "eukaryota":
        table_name = 'busco_eukaryota_database'
    else:
        table_name = 'busco_virus_database'
    add_data_table(data_tables, table_name)

    # Fetch data from specified data sources
    print("Download and build database")
    download_busco_db(
        data_tables,
        table_name,
        args.database,
        args.index,
        args.version,
        target_dp)

    # Write output JSON
    print("Outputting JSON")
    with open(json_fp, 'w') as fh:
        json.dump(data_tables, fh, sort_keys=True)
    print("Done.")
