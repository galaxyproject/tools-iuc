#!/usr/bin/env python
#
# Data manager for reference data for the 'humann2' Galaxy tools
import datetime
import json
import optparse
import os
import shutil
import subprocess
import sys


HUMANN2_REFERENCE_DATA = {
    "full": "Full",
    "DEMO": "Demo",
    "uniref50_diamond": "Full UniRef50",
    "uniref50_ec_filtered_diamond": "EC-filtered UniRef50",
    "uniref50_GO_filtered_rapsearch2": "GO filtered UniRef50 for rapsearch2",
    "uniref90_diamond": "Full UniRef90",
    "uniref90_ec_filtered_diamond": "EC-filtered UniRef90",
    "DEMO_diamond": "Demo"
}


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


def download_humann2_db(data_tables, table_name, database, build, target_dir):
    """Download HUMAnN2 database

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      table_name: name of the table
      database: database to download (chocophlan or uniref)
      build: build of the database to download
      target_dir: directory to put copy or link to the data file

    """
    value = "%s-%s-%s" % (database, build, datetime.date.today().isoformat())
    db_target_dir = os.path.join(target_dir, database)
    build_target_dir = os.path.join(db_target_dir, build)
    cmd = "humann2_databases --download %s %s %s --update-config no" % (
        database,
        build,
        db_target_dir)
    subprocess.check_call(cmd, shell=True)
    shutil.move(os.path.join(db_target_dir, database), build_target_dir)
    add_data_table_entry(
        data_tables,
        table_name,
        dict(
            dbkey=build,
            value=value,
            name=HUMANN2_REFERENCE_DATA[build],
            path=build_target_dir))


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = optparse.OptionParser(description='Download HUMAnN2 database')
    parser.add_option('--database', help="Database name")
    parser.add_option('--build', help="Build of the database")
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
    os.mkdir(target_dir, mode=0o755)

    # Set up data tables dictionary
    data_tables = create_data_tables_dict()

    if options.database == "chocophlan":
        table_name = 'humann2_nucleotide_database'
    else:
        table_name = 'humann2_protein_database'
    add_data_table(data_tables, table_name)

    # Fetch data from specified data sources
    download_humann2_db(
        data_tables,
        table_name,
        options.database,
        options.build,
        target_dir)

    # Write output JSON
    print("Outputting JSON")
    print(str(json.dumps(data_tables)))
    open(jsonfile, 'wb').write(json.dumps(data_tables))
    print("Done.")
