#!/usr/bin/env python
#
# Data manager for reference data for the 'mothur_toolsuite' Galaxy tools
import json
import optparse
import os
import shutil
import sys
import tarfile
import tempfile
import urllib2
import zipfile

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
        "GS FLX Titanium": ["http://www.mothur.org/w/images/9/96/LookUp_Titanium.zip", ]
    },
    "lookup_gsflx": {
        "GSFLX": ["http://www.mothur.org/w/images/8/84/LookUp_GSFLX.zip", ]
    },
    "lookup_gs20": {
        "GS20": ["http://www.mothur.org/w/images/7/7b/LookUp_GS20.zip", ]
    },
    # RDP reference files
    # http://www.mothur.org/wiki/RDP_reference_files
    "RDP_v16": {
        "16S rRNA RDP training set 16":
        ["https://mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz", ],
        "16S rRNA PDS training set 16":
        ["https://mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz", ],
    },
    "RDP_v14": {
        "16S rRNA RDP training set 14":
        ["https://mothur.org/w/images/6/6c/Trainset14_032015.rdp.tgz", ],
        "16S rRNA PDS training set 14":
        ["https://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz", ],
    },
    "RDP_v10": {
        "16S rRNA RDP training set 10":
        ["http://www.mothur.org/w/images/b/b5/Trainset10_082014.rdp.tgz", ],
        "16S rRNA PDS training set 10":
        ["http://www.mothur.org/w/images/2/24/Trainset10_082014.pds.tgz", ],
    },
    "RDP_v9": {
        "16S rRNA RDP training set 9":
        ["http://www.mothur.org/w/images/7/72/Trainset9_032012.rdp.zip", ],
        "16S rRNA PDS training set 9":
        ["http://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip", ],
    },
    "RDP_v7": {
        "16S rRNA RDP training set 7":
        ["http://www.mothur.org/w/images/2/29/Trainset7_112011.rdp.zip", ],
        "16S rRNA PDS training set 7":
        ["http://www.mothur.org/w/images/4/4a/Trainset7_112011.pds.zip", ],
        "8S rRNA Fungi training set 7":
        ["http://www.mothur.org/w/images/3/36/FungiLSU_train_v7.zip", ],
    },
    "RDP_v6": {
        "RDP training set 6":
        ["http://www.mothur.org/w/images/4/49/RDPTrainingSet.zip", ],
    },
    # Silva reference files
    # http://www.mothur.org/wiki/Silva_reference_files
    "silva_release_128": {
        "SILVA release 128":
        ["https://mothur.org/w/images/b/b4/Silva.nr_v128.tgz",
         "https://mothur.org/w/images/a/a4/Silva.seed_v128.tgz", ],
    },
    "silva_release_123": {
        "SILVA release 123":
        ["https://mothur.org/w/images/b/be/Silva.nr_v123.tgz",
         "https://mothur.org/w/images/1/15/Silva.seed_v123.tgz", ],
    },
    "silva_release_119": {
        "SILVA release 119":
        ["http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz",
         "http://www.mothur.org/w/images/5/56/Silva.seed_v119.tgz", ],
    },
    "silva_release_102": {
        "SILVA release 102":
        ["http://www.mothur.org/w/images/9/98/Silva.bacteria.zip",
         "http://www.mothur.org/w/images/3/3c/Silva.archaea.zip",
         "http://www.mothur.org/w/images/1/1a/Silva.eukarya.zip", ],
    },
    "silva_gold_bacteria": {
        "SILVA gold":
        ["http://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip", ],
    },
    # Greengenes
    # http://www.mothur.org/wiki/Greengenes-formatted_databases
    "greengenes_August2013": {
        "Greengenes August 2013":
        ["http://www.mothur.org/w/images/1/19/Gg_13_8_99.refalign.tgz",
         "http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz", ],
    },
    "greengenes_May2013": {
        "Greengenes May 2013":
        ["http://www.mothur.org/w/images/c/cd/Gg_13_5_99.refalign.tgz",
         "http://www.mothur.org/w/images/9/9d/Gg_13_5_99.taxonomy.tgz", ],
    },
    "greengenes_old": {
        "Greengenes pre-May 2013":
        ["http://www.mothur.org/w/images/7/72/Greengenes.alignment.zip",
         "http://www.mothur.org/w/images/1/16/Greengenes.tax.tgz", ],
    },
    "greengenes_gold_alignment": {
        "Greengenes gold alignment":
        ["http://www.mothur.org/w/images/2/21/Greengenes.gold.alignment.zip", ],
    },
    # Secondary structure maps
    # http://www.mothur.org/wiki/Secondary_structure_map
    "secondary_structure_maps_silva": {
        "SILVA":
        ["http://www.mothur.org/w/images/6/6d/Silva_ss_map.zip", ],
    },
    "secondary_structure_maps_greengenes": {
        "Greengenes":
        ["http://www.mothur.org/w/images/4/4b/Gg_ss_map.zip", ],
    },
    # Lane masks: not used here?
    "lane_masks": {
        "Greengenes-compatible":
        ["http://www.mothur.org/w/images/2/2a/Lane1241.gg.filter",
         "http://www.mothur.org/w/images/a/a0/Lane1287.gg.filter",
         "http://www.mothur.org/w/images/3/3d/Lane1349.gg.filter", ],
        "SILVA-compatible":
        ["http://www.mothur.org/w/images/6/6d/Lane1349.silva.filter", ]
    },
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
    print "Downloading %s" % url
    if not target:
        target = os.path.basename(url)
    if wd:
        target = os.path.join(wd, target)
    print "Saving to %s" % target
    open(target, 'wb').write(urllib2.urlopen(url).read())
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
        print "%s: not ZIP formatted file"
        return [filen]
    file_list = []
    z = zipfile.ZipFile(filen)
    for name in z.namelist():
        if reduce(lambda x, y: x or name.startswith(y), IGNORE_PATHS, False):
            print "Ignoring %s" % name
            continue
        if wd:
            target = os.path.join(wd, name)
        else:
            target = name
        if name.endswith('/'):
            # Make directory
            print "Creating dir %s" % target
            try:
                os.makedirs(target)
            except OSError:
                pass
        else:
            # Extract file
            print "Extracting %s" % name
            try:
                os.makedirs(os.path.dirname(target))
            except OSError:
                pass
            open(target, 'wb').write(z.read(name))
            file_list.append(target)
    print "Removing %s" % filen
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
        print "%s: not TAR file"
        return [filen]
    t = tarfile.open(filen)
    for name in t.getnames():
        # Check for unwanted files
        if reduce(lambda x, y: x or name.startswith(y), IGNORE_PATHS, False):
            print "Ignoring %s" % name
            continue
        # Extract file
        print "Extracting %s" % name
        t.extract(name, wd)
        if wd:
            target = os.path.join(wd, name)
        else:
            target = name
        file_list.append(target)
    print "Removing %s" % filen
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
    print "Unpack %s" % filen
    ext = os.path.splitext(filen)[1]
    print "Extension: %s" % ext
    if ext == ".zip":
        return unpack_zip_archive(filen, wd=wd)
    elif ext == ".tgz":
        return unpack_tar_archive(filen, wd=wd)
    else:
        return [filen]


def fetch_files(urls, wd=None, files=None):
    """Download and unpack files from a list of URLs

    Given a list of URLs, download and unpack each
    one, and return a list of the extracted files.

    'wd' specifies the working directory to extract
    the files to, otherwise they are extracted to the
    current working directory.

    If 'files' is given then the list of extracted
    files will be appended to this list before being
    returned.

    """
    if files is None:
        files = []
    for url in urls:
        filen = download_file(url, wd=wd)
        files.extend(unpack_archive(filen, wd=wd))
    return files


# Utility functions specific to the Mothur reference data
def identify_type(filen):
    """Return the data table name based on the file name

    """
    ext = os.path.splitext(filen)[1]
    try:
        return MOTHUR_FILE_TYPES[ext]
    except KeyError:
        print "WARNING: unknown file type for " + filen + ", skipping"
        return None


def get_name(filen):
    """Generate a descriptive name based on the file name
    """
    # type_ = identify_type(filen)
    name = os.path.splitext(os.path.basename(filen))[0]
    for delim in ('.', '_'):
        name = name.replace(delim, ' ')
    return name


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
    print "Working dir %s" % wd
    # Iterate over all requested reference data URLs
    for dataset in datasets:
        print "Handling dataset '%s'" % dataset
        for name in MOTHUR_REFERENCE_DATA[dataset]:
            for f in fetch_files(MOTHUR_REFERENCE_DATA[dataset][name], wd=wd):
                type_ = identify_type(f)
                entry_name = "%s (%s)" % (os.path.splitext(os.path.basename(f))[0], name)
                print "%s\t\'%s'\t.../%s" % (type_, entry_name, os.path.basename(f))
                if type_ is not None:
                    # Move to target dir
                    ref_data_file = os.path.basename(f)
                    f1 = os.path.join(target_dir, ref_data_file)
                    print "Moving %s to %s" % (f, f1)
                    os.rename(f, f1)
                    # Add entry to data table
                    table_name = "mothur_%s" % type_
                    add_data_table_entry(data_tables, table_name, dict(name=entry_name, value=ref_data_file))
    # Remove working dir
    print "Removing %s" % wd
    shutil.rmtree(wd)


def files_from_filesystem_paths(paths):
    """Return list of file paths from arbitrary input paths

    Given a list of filesystem paths, return a list of
    full paths corresponding to all files found recursively
    from under those paths.

    """
    # Collect files to add
    files = []
    for path in paths:
        path = os.path.abspath(path)
        print "Examining '%s'..." % path
        if os.path.isfile(path):
            # Store full path for file
            files.append(path)
        elif os.path.isdir(path):
            # Descend into directory and collect the files
            for f in os.listdir(path):
                files.extend(files_from_filesystem_paths((os.path.join(path, f), )))
        else:
            print "Not a file or directory, ignored"
    return files


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
    # Collect list of files based on input paths
    files = files_from_filesystem_paths(paths)
    # Handle each file individually
    for f in files:
        type_ = identify_type(f)
        if type_ is None:
            print "%s: unrecognised type, skipped" % f
            continue
        ref_data_file = os.path.basename(f)
        target_file = os.path.join(target_dir, ref_data_file)
        entry_name = "%s" % os.path.splitext(ref_data_file)[0]
        if description:
            entry_name += " (%s)" % description
        print "%s\t\'%s'\t.../%s" % (type_, entry_name, ref_data_file)
        # Link to or copy the data
        if link_to_data:
            os.symlink(f, target_file)
        else:
            shutil.copyfile(f, target_file)
        # Add entry to data table
        table_name = "mothur_%s" % type_
        add_data_table_entry(data_tables, table_name, dict(name=entry_name, value=ref_data_file))


if __name__ == "__main__":
    print "Starting..."

    # Read command line
    parser = optparse.OptionParser()
    parser.add_option('--source', action='store', dest='data_source')
    parser.add_option('--datasets', action='store', dest='datasets', default='')
    parser.add_option('--paths', action='store', dest='paths', default=[])
    parser.add_option('--description', action='store', dest='description', default='')
    parser.add_option('--link', action='store_true', dest='link_to_data')
    options, args = parser.parse_args()
    print "options: %s" % options
    print "args   : %s" % args

    # Check for JSON file
    if len(args) != 1:
        sys.stderr.write("Need to supply JSON file name")
        sys.exit(1)

    jsonfile = args[0]

    # Read the input JSON
    params, target_dir = read_input_json(jsonfile)

    # Make the target directory
    print "Making %s" % target_dir
    os.mkdir(target_dir)

    # Set up data tables dictionary
    data_tables = create_data_tables_dict()
    add_data_table(data_tables, 'mothur_lookup')
    add_data_table(data_tables, 'mothur_aligndb')
    add_data_table(data_tables, 'mothur_map')
    add_data_table(data_tables, 'mothur_taxonomy')

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
    print "Outputting JSON"
    print str(json.dumps(data_tables))
    open(jsonfile, 'wb').write(json.dumps(data_tables))
    print "Done."
