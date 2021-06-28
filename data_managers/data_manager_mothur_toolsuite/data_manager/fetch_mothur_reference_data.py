#!/usr/bin/env python3
#
# Data manager for reference data for the 'mothur_toolsuite' Galaxy tools
import argparse
import io
import json
import os
import pathlib
import shutil
import tarfile
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from functools import reduce

# When extracting files from archives, skip names that
# start with the following strings
IGNORE_PATHS = ('.', '__MACOSX/', '__')

# Map file extensions to data table names
MOTHUR_FILE_TYPES = {".map": "map",
                     ".fasta": "aligndb",
                     ".align": "aligndb",
                     ".pat": "lookup",
                     ".tax": "taxonomy"}

# Reference data URLs
MOTHUR_REFERENCE_DATA = {
    # Look up data
    # http://www.mothur.org/wiki/Lookup_files
    "lookup_titanium": {
        "GS FLX Titanium": ["https://mothur.s3.us-east-2.amazonaws.com/wiki/lookup_titanium.zip", ]
    },
    "lookup_gsflx": {
        "GSFLX": ["https://mothur.s3.us-east-2.amazonaws.com/wiki/lookup_gsflx.zip", ]
    },
    "lookup_gs20": {
        "GS20": ["https://mothur.s3.us-east-2.amazonaws.com/wiki/lookup_gs20.zip", ]
    },
    # RDP reference files
    # http://www.mothur.org/wiki/RDP_reference_files
    "RDP_v18": {
        "16S rRNA RDP training set 18":
            [
                "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.rdp.tgz", ],
        "16S rRNA PDS training set 18":
            [
                "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.pds.tgz", ],
    },
    "RDP_v16": {
        "16S rRNA RDP training set 16":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.rdp.tgz", ],
        "16S rRNA PDS training set 16":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz", ],
    },
    "RDP_v14": {
        "16S rRNA RDP training set 14":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset14_032015.rdp.tgz", ],
        "16S rRNA PDS training set 14":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset14_032015.pds.tgz", ],
    },
    "RDP_v10": {
        "16S rRNA RDP training set 10":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset10_082014.rdp.tgz", ],
        "16S rRNA PDS training set 10":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset10_082014.pds.tgz", ],
    },
    "RDP_v9": {
        "16S rRNA RDP training set 9":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.rdp.zip", ],
        "16S rRNA PDS training set 9":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip", ],
    },
    "RDP_v7": {
        "16S rRNA RDP training set 7":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset7_112011.rdp.zip", ],
        "16S rRNA PDS training set 7":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset7_112011.pds.zip", ],
        "8S rRNA Fungi training set 7":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/fungilsu_train_v7.zip", ],
    },
    "RDP_v6": {
        "RDP training set 6":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/rdptrainingset.zip", ],
    },
    # Silva reference files
    # http://www.mothur.org/wiki/Silva_reference_files
    "silva_release_138.1": {
        "SILVA release 138.1":
            [
                "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz",
                "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz", ],
    },
    "silva_release_128": {
        "SILVA release 128":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v128.tgz",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v128.tgz", ],
    },
    "silva_release_123": {
        "SILVA release 123":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v123.tgz",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v123.tgz", ],
    },
    "silva_release_119": {
        "SILVA release 119":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v119.tgz",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v119.tgz", ],
    },
    "silva_release_102": {
        "SILVA release 102":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.archaea.zip",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.eukarya.zip", ],
    },
    "silva_gold_bacteria": {
        "SILVA gold":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip", ],
    },
    # Greengenes
    # http://www.mothur.org/wiki/Greengenes-formatted_databases
    "greengenes_August2013": {
        "Greengenes August 2013":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_8_99.refalign.tgz",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_8_99.taxonomy.tgz", ],
    },
    "greengenes_May2013": {
        "Greengenes May 2013":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_5_99.refalign.tgz",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_5_99.taxonomy.tgz", ],
    },
    "greengenes_old": {
        "Greengenes pre-May 2013":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/greengenes.alignment.zip",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/greengenes.tax.tgz", ],
    },
    "greengenes_gold_alignment": {
        "Greengenes gold alignment":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/greengenes.gold.alignment.zip", ],
    },
    # Secondary structure maps
    # http://www.mothur.org/wiki/Secondary_structure_map
    "secondary_structure_maps_silva": {
        "SILVA":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/silva_ss_map.zip", ],
    },
    "secondary_structure_maps_greengenes": {
        "Greengenes":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_ss_map.zip", ],
    },
    # Lane masks: not used here?
    "lane_masks": {
        "Greengenes-compatible":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/Lane1241.gg.filter",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/lane1287.gg.filter",
         "https://mothur.s3.us-east-2.amazonaws.com/wiki/lane1349.gg.filter", ],
        "SILVA-compatible":
        ["https://mothur.s3.us-east-2.amazonaws.com/wiki/lane1349.silva.filter", ]
    },
}


# Utility functions for downloading and unpacking archive files
def download_file(url, target=None, wd=None):
    """Download a file from a URL

    Fetches a file from the specified URL.

    If 'target' is specified then the file is saved to this
    name; otherwise it's saved as the basename of the URL.

    If 'wd' is specified then it is used as the 'working
    directory' where the file will be save on the local
    system.

    Returns the name that the file is saved with.

    """
    print(f"Downloading {url}")
    if not target:
        target = os.path.basename(url)
    if wd:
        target = os.path.join(wd, target)
    print(f"Saving to {target}")
    with open(target, 'wb') as fh:
        url_h = urllib.request.urlopen(url)
        while True:
            buffer = url_h.read(io.DEFAULT_BUFFER_SIZE)
            if buffer == b"":
                break
            fh.write(buffer)
    return target


def unpack_zip_archive(filen, wd=None):
    """Extract files from a ZIP archive

    Given a ZIP archive, extract the files it contains
    and return a list of the resulting file names and
    paths.

    'wd' specifies the working directory to extract
    the files to, otherwise they are extracted to the
    current working directory.

    Once all the files are extracted the ZIP archive
    file is deleted from the file system.

    """
    if not zipfile.is_zipfile(filen):
        print(f"{filen}: not ZIP formatted file")
        return [filen]
    file_list = []
    with zipfile.ZipFile(filen) as z:
        for name in z.namelist():
            if reduce(lambda x, y: x or name.startswith(y), IGNORE_PATHS, False):
                print(f"Ignoring {name}")
                continue
            target = os.path.join(wd, name) if wd else name
            if name.endswith('/'):
                # Make directory
                print(f"Creating dir {target}")
                os.makedirs(target, exist_ok=True)
            else:
                # Extract file
                print("Extracting {target}")
                os.makedirs(os.path.dirname(target), exist_ok=True)
                with open(target, 'wb') as fh:
                    fh.write(z.read(name))
                file_list.append(target)
    print(f"Removing {filen}")
    os.remove(filen)
    return file_list


def unpack_tar_archive(filen, wd=None):
    """Extract files from a TAR archive

    Given a TAR archive (which optionally can be
    compressed with either gzip or bz2), extract the
    files it contains and return a list of the
    resulting file names and paths.

    'wd' specifies the working directory to extract
    the files to, otherwise they are extracted to the
    current working directory.

    Once all the files are extracted the TAR archive
    file is deleted from the file system.

    """
    file_list = []
    if not tarfile.is_tarfile(filen):
        print(f"{filen}: not TAR file")
        return [filen]
    with tarfile.open(filen) as t:
        for name in t.getnames():
            # Check for unwanted files
            if reduce(lambda x, y: x or name.startswith(y), IGNORE_PATHS, False):
                print(f"Ignoring {name}")
                continue
            # Extract file
            print(f"Extracting {name}")
            t.extract(name, wd)
            target = os.path.join(wd, name) if wd else name
            file_list.append(target)
    print(f"Removing {filen}")
    os.remove(filen)
    return file_list


def unpack_archive(filen, wd=None):
    """Extract files from an archive

    Wrapper function that calls the appropriate
    unpacking function depending on the archive
    type, and returns a list of files that have
    been extracted.

    'wd' specifies the working directory to extract
    the files to, otherwise they are extracted to the
    current working directory.

    """
    print(f"Unpack {filen}")
    ext = os.path.splitext(filen)[1]
    print(f"Extension: {ext}")
    if ext == ".zip":
        return unpack_zip_archive(filen, wd=wd)
    elif ext == ".tgz":
        return unpack_tar_archive(filen, wd=wd)
    else:
        return [filen]


# Utility functions specific to the Mothur reference data
def identify_type(filen):
    """Return the data table name based on the file name

    """
    ext = os.path.splitext(filen)[1]
    try:
        return MOTHUR_FILE_TYPES[ext]
    except KeyError:
        print(f"WARNING: unknown file type for {filen}, skipping")
        return None


def fetch_from_mothur_website(data_tables, target_dir, datasets):
    """Fetch reference data from the Mothur website

    For each dataset in the list 'datasets', download (and if
    necessary unpack) the related files from the Mothur website,
    copy them to the data manager's target directory, and add
    references to the files to the appropriate data table.

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      target_dir: directory to put the downloaded files
      datasets: a list of dataset names corresponding to keys in
        the MOTHUR_REFERENCE_DATA dictionary
    """
    # Make working dir
    wd = tempfile.mkdtemp(suffix=".mothur", dir=os.getcwd())
    print(f"Working dir {wd}")
    # Iterate over all requested reference data URLs
    for dataset in datasets:
        print(f"Handling dataset '{dataset}'")
        for name, urls in MOTHUR_REFERENCE_DATA[dataset].items():
            for url in urls:
                filen = download_file(url, wd=wd)
                for f in unpack_archive(filen, wd=wd):
                    type_ = identify_type(f)
                    name_from_file = os.path.splitext(os.path.basename(f))[0]
                    entry_name = f"{name_from_file} ({name})"
                    print(f"{type_}\t\'{entry_name}'\t.../{os.path.basename(f)}")
                    if type_ is not None:
                        # Move to target dir
                        ref_data_file = os.path.basename(f)
                        f1 = os.path.join(target_dir, ref_data_file)
                        print(f"Moving {f} to {f1}")
                        shutil.move(f, f1)
                        data_tables['data_tables'][f"mothur_{type_}"].append(
                            dict(name=entry_name, value=ref_data_file))
    # Remove working dir
    print(f"Removing {wd}")
    shutil.rmtree(wd)


def files_from_filesystem_paths(paths):
    """Return list of file paths from arbitrary input paths

    Given a list of filesystem paths, return a list of
    full paths corresponding to all files found recursively
    from under those paths.

    """
    # Collect files to add
    for raw_path in paths:
        path = pathlib.Path(raw_path).absolute()
        print(f"Examining '{str(path)}'...")
        if path.is_file():
            yield str(path)
        elif path.is_dir():
            yield from files_from_filesystem_paths(path.iterdir())
        else:
            print("Not a file or directory, ignored")


def import_from_server(data_tables, target_dir, paths, description, link_to_data=False):
    """Import reference data from filesystem paths

    Creates references to the specified file(s) on the Galaxy
    server in the appropriate data table (determined from the
    file extension).

    The 'data_tables' dictionary should have been created using
    the 'create_data_tables_dict' and 'add_data_table' functions.

    Arguments:
      data_tables: a dictionary containing the data table info
      target_dir: directory to put copy or link to the data file
      paths: list of file and/or directory paths to import
      description: text to associate with the files
      link_to_data: boolean, if False then copy the data file
        into Galaxy (default); if True then make a symlink to
        the data file

    """
    for f in files_from_filesystem_paths(paths):
        type_ = identify_type(f)
        if type_ is None:
            print(f"{f}: unrecognised type, skipped")
            continue
        ref_data_file = os.path.basename(f)
        target_file = os.path.join(target_dir, ref_data_file)
        entry_name = os.path.splitext(ref_data_file)[0]
        if description:
            entry_name += f" ({description})"
        print(f"{type_}\t\'{entry_name}'\t.../{ref_data_file}")
        # Link to or copy the data
        if link_to_data:
            os.symlink(f, target_file)
        else:
            shutil.copyfile(f, target_file)
        data_tables['data_tables'][f"mothur_{type_}"].append(
            dict(name=entry_name, value=ref_data_file))


if __name__ == "__main__":
    print("Starting...")

    # Read command line
    parser = argparse.ArgumentParser()
    parser.add_argument("json", metavar="JSON")
    parser.add_argument('--source', action='store', dest='data_source')
    parser.add_argument('--datasets', action='store', dest='datasets', default='')
    parser.add_argument('--paths', action='store', dest='paths', default=[])
    parser.add_argument('--description', action='store', dest='description', default='')
    parser.add_argument('--link', action='store_true', dest='link_to_data')
    options = parser.parse_args()

    with open(options.json, "rt") as fh:
        json_params = json.load(fh)
    params = json_params["param_dict"]
    target_dir = json_params["output_data"][0]["extra_files_path"]

    # Make the target directory
    print(f"Making {target_dir}")
    os.mkdir(target_dir)

    # Set up data tables dictionary
    data_tables = {
        "data_tables": {
            "mothur_lookup": [],
            "mothur_aligndb": [],
            "mothur_map": [],
            "mothur_taxonomy": []
        }}

    # Fetch data from specified data sources
    if options.data_source == 'mothur_website':
        datasets = options.datasets.split(',')
        fetch_from_mothur_website(data_tables, target_dir, datasets)
    elif options.data_source == 'filesystem_paths':
        # Check description text
        description = options.description.strip()
        # Get list of paths (need to remove any escapes for '\n' and '\r'
        # that might have been inserted by Galaxy)
        paths = options.paths.replace('__cn__', '\n').replace('__cr__', '\r').split()
        import_from_server(data_tables, target_dir, paths, description, link_to_data=options.link_to_data)
    # Write output JSON
    print("Outputting JSON")
    with open(options.json, 'w') as fh:
        json.dump(data_tables, fh, sort_keys=True)
    print("Done.")
