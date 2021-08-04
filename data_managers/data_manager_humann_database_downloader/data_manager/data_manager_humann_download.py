#!/usr/bin/env python
#
# Data manager for reference data for the 'humann' Galaxy tools
import argparse
import json
import subprocess
from datetime import date
from pathlib import Path

HUMANN_REFERENCE_DATA = {
    "chocophlan": {
        "full": "Full ChocoPhlAn for HUManN",
        "DEMO": "Demo ChocoPhlAn for HUManN"
    },
    "uniref": {
        "uniref50_diamond": "Full UniRef50 for HUManN",
        "uniref50_ec_filtered_diamond": "EC-filtered UniRef50 for HUManN",
        "uniref90_diamond": "Full UniRef90 for HUManN",
        "uniref90_ec_filtered_diamond": "EC-filtered UniRef90 for HUManN",
        "DEMO_diamond": "Demo UniRef for HUManN"
    },
    "utility_mapping": {
        "full": {
            "map_uniref50_uniref90": "Mapping (full) for UniRef50 from UniRef90",
            "map_ko_uniref90": "Mapping (full) for KEGG Orthogroups (KOs) from UniRef90",
            "map_eggnog_name": "Mapping (full) between EggNOG (including COGs) ids and names",
            "map_uniref90_name": "Mapping (full) between UniRef90 ids and names",
            "map_go_uniref90": "Mapping (full) for Gene Ontology (GO) from UniRef90",
            "uniref90-tol-lca": "Mapping (full) for LCA for UniRef90",
            "uniref50-tol-lca": "Mapping (full) for LCA for UniRef50",
            "map_eggnog_uniref50": "Mapping (full) for EggNOG (including COGs) from UniRef50",
            "map_pfam_uniref90": "Mapping (full) for Pfam domains from UniRef90",
            "map_go_uniref50": "Mapping (full) for Gene Ontology (GO) from UniRef50",
            "map_ko_name": "Mapping (full) between KEGG Orthogroups (KOs) ids and names",
            "map_level4ec_uniref90": "Mapping (full) for Level-4 enzyme commission (EC) categories from UniRef90",
            "map_go_name": "Mapping (full) between Gene Ontology (GO) ids and names",
            "map_ko_uniref50": "Mapping (full) for KEGG Orthogroups (KOs) from UniRef50",
            "map_level4ec_uniref50": "Mapping (full) for Level-4 enzyme commission (EC) categories from UniRef90",
            "map_pfam_uniref50": "Mapping (full) for Pfam domains from UniRef50",
            "map_eggnog_uniref90": "Mapping (full) for EggNOG (including COGs) from UniRef90",
            "map_uniref50_name": "Mapping (full) between UniRef50 ids and names",
            "map_ec_name": "Mapping (full) between Level-4 enzyme commission (EC) categories ids and names",
            "map_pfam_name": "Mapping (full) between Pfam domains ids and names"
        }
    }
}


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


def download_humann_db(data_tables, table_name, database, build, version, target_dp):
    """Download HUMAnN database

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
      version: tool version
      target_dp: directory to put copy or link to the data file
    """
    db_target_dp = target_dp / Path(database)
    db_dp = db_target_dp / Path(database)
    build_target_dp = db_target_dp / Path(build)
    # launch tool to get db
    cmd = "humann_databases --download %s %s %s --update-config no" % (
        database,
        build,
        db_target_dp)
    subprocess.check_call(cmd, shell=True)
    # move db
    db_dp.rename(build_target_dp)
    # add details to data table
    if database != "utility_mapping":
        add_data_table_entry(
            data_tables,
            table_name,
            dict(
                value="%s-%s-%s-%s" % (database, build, version, date.today().strftime("%d%m%Y")),
                name=HUMANN_REFERENCE_DATA[database][build],
                dbkey=version,
                path=str(build_target_dp)))
    elif args.database == "utility_mapping":
        for x in build_target_dp.iterdir():
            name = str(x.stem).split('.')[0]
            add_data_table_entry(
                data_tables,
                table_name,
                dict(
                    value="%s-%s-%s-%s-%s%s" % (database, build, name, version, date.today().strftime("%d%m%Y"), x.suffix),
                    name=HUMANN_REFERENCE_DATA["utility_mapping"][build][name],
                    dbkey=version,
                    path=str(x)))


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser(description='Download HUMAnN database')
    parser.add_argument('--database', help="Database name")
    parser.add_argument('--build', help="Build of the database")
    parser.add_argument('--version', help="HUMAnN version")
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
    if args.database == "chocophlan":
        table_name = 'humann_nucleotide_database'
    elif args.database == "uniref":
        table_name = 'humann_protein_database'
    elif args.database == "utility_mapping":
        table_name = 'humann_utility_mapping'
    add_data_table(data_tables, table_name)

    # Fetch data from specified data sources
    print("Download and build database")
    download_humann_db(
        data_tables,
        table_name,
        args.database,
        args.build,
        args.version,
        target_dp)

    # Write output JSON
    print("Outputting JSON")
    with open(json_fp, 'w') as fh:
        json.dump(data_tables, fh, sort_keys=True)
    print("Done.")
