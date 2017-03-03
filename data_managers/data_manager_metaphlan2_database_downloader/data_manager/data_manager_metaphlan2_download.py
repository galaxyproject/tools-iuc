#!/usr/bin/env python
#
# Data manager for reference data for the 'humann2' Galaxy tools
import json
import optparse
import os
import subprocess
import sys


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


def download_metaphlan2_db(data_tables, build, table_name, target_dir):
    """Download MetaPhlAn2 database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      table_name: name of the table
      target_dir: directory to put copy or link to the data file

    """
    cmd = "download_metaphlan2_db.py --output %s" % (target_dir)
    db_dir = os.path.join(target_dir, build)
    subprocess.check_call(cmd, shell=True)
    os.rename(os.path.join(target_dir, "db_v20"), db_dir)
    add_data_table_entry(
        data_tables,
        table_name,
        dict(
            dbkey=build,
            value="mpa_v20_m200",
            name="MetaPhlAn2 clade-specific marker genes",
            path=db_dir))


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = optparse.OptionParser(description='Download MetaPhlan2 database')
    parser.add_option('--database', help="Database name")
    options, args = parser.parse_args()
    print("args   : %s" % args)

    # Check for JSON file
    if len(args) != 1:
        sys.stderr.write("Need to supply JSON file name")
        sys.exit(1)

    jsonfile = args[0]

    # Read the input JSON
    params, target_dir = read_input_json(jsonfile)

    # Make the target directory
    print("Making %s" % target_dir)
    os.mkdir(target_dir)

    # Set up data tables dictionary
    data_tables = create_data_tables_dict()
    add_data_table(data_tables, "metaphlan2_database")

    # Fetch data from specified data sources
    if options.database == "db_v20":
        download_metaphlan2_db(
            data_tables,
            "v20",
            "metaphlan2_database",
            target_dir)

    # Write output JSON
    print("Outputting JSON")
    print(str(json.dumps(data_tables)))
    with open(jsonfile, 'wb') as out:
        out.write(json.dumps(data_tables))
    print("Done.")
