#!/usr/bin/env python
# Dan Blankenberg

import bz2
import gzip
import optparse
import os
import shutil
import sys
import tarfile
import tempfile
import zipfile
from ftplib import FTP
from json import dumps, loads

try:
    # For Python 3.0 and later
    from io import BytesIO as StringIO
    from io import UnsupportedOperation
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2 imports
    from StringIO import StringIO
    from urllib2 import urlopen

    UnsupportedOperation = AttributeError

CHUNK_SIZE = 2 ** 20  # 1mb


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def get_dbkey_dbname_id_name(params, dbkey_description=None):
    dbkey = params['param_dict']['dbkey_source']['dbkey']
    # TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey  # uuid.uuid4() generate and use an uuid instead?

    if params['param_dict']['dbkey_source']['dbkey_source_selector'] == 'new':
        dbkey_name = params['param_dict']['dbkey_source']['dbkey_name']
        if not dbkey_name:
            dbkey_name = dbkey
    else:
        dbkey_name = None

    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = dbkey_description
        if not sequence_name:
            sequence_name = dbkey
    return dbkey, dbkey_name, sequence_id, sequence_name


def _get_files_in_ftp_path(ftp, path):
    path_contents = []
    ftp.retrlines('MLSD %s' % (path), path_contents.append)
    return [line.split(';')[-1].lstrip() for line in path_contents]


def _get_stream_readers_for_tar(fh, tmp_dir):
    fasta_tar = tarfile.open(fileobj=fh, mode='r:*')
    return [x for x in [fasta_tar.extractfile(member) for member in fasta_tar.getmembers()] if x]


def _get_stream_readers_for_zip(fh, tmp_dir):
    """
    Unpacks all archived files in a zip file.
    Individual files will be concatenated (in _stream_fasta_to_file)
    """
    fasta_zip = zipfile.ZipFile(fh, 'r')
    rval = []
    for member in fasta_zip.namelist():
        fasta_zip.extract(member, tmp_dir)
        rval.append(open(os.path.join(tmp_dir, member), 'rb'))
    return rval


def _get_stream_readers_for_gzip(fh, tmp_dir):
    return [gzip.GzipFile(fileobj=fh, mode='rb')]


def _get_stream_readers_for_bz2(fh, tmp_dir):
    return [bz2.BZ2File(fh.name, 'rb')]


def sort_fasta(fasta_filename, sort_method, params):
    if sort_method is None:
        return
    assert sort_method in SORTING_METHODS, ValueError("%s is not a valid sorting option." % sort_method)
    return SORTING_METHODS[sort_method](fasta_filename, params)


def _move_and_index_fasta_for_sorting(fasta_filename):
    unsorted_filename = tempfile.NamedTemporaryFile().name
    shutil.move(fasta_filename, unsorted_filename)
    fasta_offsets = {}
    unsorted_fh = open(unsorted_filename)
    while True:
        offset = unsorted_fh.tell()
        line = unsorted_fh.readline()
        if not line:
            break
        if line.startswith(">"):
            line = line.split(None, 1)[0][1:]
            fasta_offsets[line] = offset
    unsorted_fh.close()
    current_order = [_[1] for _ in sorted((_[1], _[0]) for _ in fasta_offsets.items())]
    return (unsorted_filename, fasta_offsets, current_order)


def _write_sorted_fasta(sorted_names, fasta_offsets, sorted_fasta_filename, unsorted_fasta_filename):
    unsorted_fh = open(unsorted_fasta_filename)
    sorted_fh = open(sorted_fasta_filename, 'wb+')

    for name in sorted_names:
        offset = fasta_offsets[name]
        unsorted_fh.seek(offset)
        sorted_fh.write(unsorted_fh.readline())
        while True:
            line = unsorted_fh.readline()
            if not line or line.startswith(">"):
                break
            sorted_fh.write(line)
    unsorted_fh.close()
    sorted_fh.close()


def _sort_fasta_as_is(fasta_filename, params):
    return


def _sort_fasta_lexicographical(fasta_filename, params):
    (unsorted_filename, fasta_offsets, current_order) = _move_and_index_fasta_for_sorting(fasta_filename)
    sorted_names = sorted(fasta_offsets.keys())
    if sorted_names == current_order:
        shutil.move(unsorted_filename, fasta_filename)
    else:
        _write_sorted_fasta(sorted_names, fasta_offsets, fasta_filename, unsorted_filename)


def _sort_fasta_gatk(fasta_filename, params):
    (unsorted_filename, fasta_offsets, current_order) = _move_and_index_fasta_for_sorting(fasta_filename)
    sorted_names = list(map(str, range(1, 23))) + ['X', 'Y']
    # detect if we have chrN, or just N
    has_chr = False
    for chrom in sorted_names:
        if "chr%s" % chrom in current_order:
            has_chr = True
            break

    if has_chr:
        sorted_names = ["chr%s" % x for x in sorted_names]
        sorted_names.insert(0, "chrM")
    else:
        sorted_names.insert(0, "MT")
    sorted_names.extend("%s_random" % x for x in sorted_names)

    existing_sorted_names = []
    for name in sorted_names:
        if name in current_order:
            existing_sorted_names.append(name)
    for name in current_order:
        # TODO: confirm that non-canonical names do not need to be sorted specially
        if name not in existing_sorted_names:
            existing_sorted_names.append(name)

    if existing_sorted_names == current_order:
        shutil.move(unsorted_filename, fasta_filename)
    else:
        _write_sorted_fasta(existing_sorted_names, fasta_offsets, fasta_filename, unsorted_filename)


def _sort_fasta_custom(fasta_filename, params):
    (unsorted_filename, fasta_offsets, current_order) = _move_and_index_fasta_for_sorting(fasta_filename)
    sorted_names = []
    for id_repeat in params['param_dict']['sorting']['sequence_identifiers']:
        sorted_names.append(id_repeat['identifier'])
    handle_not_listed = params['param_dict']['sorting']['handle_not_listed_selector']
    if handle_not_listed.startswith('keep'):
        add_list = []
        for name in current_order:
            if name not in sorted_names:
                add_list.append(name)
        if add_list:
            if handle_not_listed == 'keep_append':
                sorted_names.extend(add_list)
            else:
                add_list.extend(sorted_names)
                sorted_names = add_list
    if sorted_names == current_order:
        shutil.move(unsorted_filename, fasta_filename)
    else:
        _write_sorted_fasta(sorted_names, fasta_offsets, fasta_filename, unsorted_filename)


def _download_file(start, fh):
    tmp = tempfile.NamedTemporaryFile()
    tmp.write(start)
    while True:
        data = fh.read(CHUNK_SIZE)
        if data:
            tmp.write(data)
        else:
            break
    tmp.flush()
    tmp.seek(0)
    return tmp


def get_stream_reader(fh, tmp_dir):
    """
    Check if file is compressed and return correct stream reader.
    If file has to be downloaded, do it now.
    """
    magic_dict = {
        b"\x1f\x8b\x08": _get_stream_readers_for_gzip,
        b"\x42\x5a\x68": _get_stream_readers_for_bz2,
        b"\x50\x4b\x03\x04": _get_stream_readers_for_zip,
    }
    start_of_file = fh.read(CHUNK_SIZE)
    try:
        fh.seek(0)
    except UnsupportedOperation:  # This happens if fh has been created by urlopen
        fh = _download_file(start_of_file, fh)
    try:  # Check if file is tar file
        if tarfile.open(fileobj=StringIO(start_of_file)):
            return _get_stream_readers_for_tar(fh, tmp_dir)
    except tarfile.ReadError:
        pass
    for k, v in magic_dict.items():
        if start_of_file.startswith(k):
            return v(fh, tmp_dir)
    return [fh]


def _get_ucsc_download_address(params, dbkey):
    """
    Check if we can find the correct file for the supplied dbkey on UCSC's FTP server
    """
    UCSC_FTP_SERVER = 'hgdownload.cse.ucsc.edu'
    UCSC_DOWNLOAD_PATH = '/goldenPath/%s/bigZips/'
    COMPRESSED_EXTENSIONS = ['.tar.gz', '.tgz', '.tar.bz2', '.zip', '.fa.gz', '.fa.bz2']

    email = params['param_dict']['__user_email__']
    if not email:
        email = 'anonymous@example.com'

    ucsc_dbkey = params['param_dict']['reference_source']['requested_dbkey'] or dbkey
    UCSC_CHROM_FA_FILENAMES = ['%s.chromFa' % ucsc_dbkey, 'chromFa', ucsc_dbkey]

    ftp = FTP(UCSC_FTP_SERVER)
    ftp.login('anonymous', email)

    ucsc_path = UCSC_DOWNLOAD_PATH % ucsc_dbkey
    path_contents = _get_files_in_ftp_path(ftp, ucsc_path)
    ftp.quit()

    for ucsc_chrom_fa_filename in UCSC_CHROM_FA_FILENAMES:
        for ext in COMPRESSED_EXTENSIONS:
            if "%s%s" % (ucsc_chrom_fa_filename, ext) in path_contents:
                ucsc_file_name = "%s%s%s" % (ucsc_path, ucsc_chrom_fa_filename, ext)
                return "ftp://%s%s" % (UCSC_FTP_SERVER, ucsc_file_name)

    raise Exception('Unable to determine filename for UCSC Genome for %s: %s' % (ucsc_dbkey, path_contents))


def add_fasta_to_table(data_manager_dict, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params, **kwds):
    fasta_filename = _stream_fasta_to_file(fasta_stream=fasta_readers,
                                           target_directory=target_directory,
                                           sequence_id=sequence_id)
    sort_fasta(fasta_filename, params['param_dict']['sorting']['sort_selector'], params)
    if dbkey_name:
        add_dbkey_to_table(data_manager_dict=data_manager_dict,
                           target_directory=target_directory,
                           dbkey=dbkey,
                           dbkey_name=dbkey_name,
                           fasta_filename=fasta_filename)
    _add_data_table_entry(data_manager_dict,
                          data_table_entry=dict(value=sequence_id, dbkey=dbkey, name=sequence_name, path=os.path.basename(fasta_filename)),
                          data_table_name='all_fasta')
    return fasta_filename


def add_dbkey_to_table(data_manager_dict, target_directory, dbkey, dbkey_name, fasta_filename):
    # do len calc here
    len_base_name = "%s.len" % (dbkey)
    compute_fasta_length(fasta_filename, os.path.join(target_directory, len_base_name), keep_first_word=True)
    dbkey_dict = dict(value=dbkey, name=dbkey_name, len_path=len_base_name)
    _add_data_table_entry(data_manager_dict,
                          data_table_entry=dbkey_dict,
                          data_table_name='__dbkeys__')


def download_from_ucsc(params, dbkey, tmp_dir, **kwds):
    url = _get_ucsc_download_address(params, dbkey)
    return get_stream_reader(urlopen(url), tmp_dir)


def download_from_ncbi(params, tmp_dir, **kwds):
    NCBI_DOWNLOAD_URL = 'http://togows.dbcls.jp/entry/ncbi-nucleotide/%s.fasta'  # FIXME: taken from dave's genome manager...why some japan site?
    requested_identifier = params['param_dict']['reference_source']['requested_identifier']
    url = NCBI_DOWNLOAD_URL % requested_identifier
    return get_stream_reader(urlopen(url), tmp_dir)


def download_from_url(params, tmp_dir, **kwds):
    """
    Download a file from a URL and return a list of filehandles from which to read the data.

    >>> url = 'https://github.com/mvdbeek/tools-devteam/raw/data_manager/data_managers/data_manager_fetch_genome_dbkeys_all_fasta/test-data/test.tar'
    >>> params = {'param_dict': {'reference_source': {'user_url': url}}}
    >>> tmp_dir = tempfile.mkdtemp()
    >>> fh = download_from_url(params=params, tmp_dir=tmp_dir)[0][0]
    >>> assert fh.readline().startswith('>FBtr0304171')
    >>> url = 'https://github.com/mvdbeek/tools-devteam/raw/data_manager/data_managers/data_manager_fetch_genome_dbkeys_all_fasta/test-data/test.tar.bz2'
    >>> params = {'param_dict': {'reference_source': {'user_url': url}}}
    >>> fh = download_from_url(params=params, tmp_dir=tmp_dir)[0][0]
    >>> assert fh.readline().startswith('>FBtr0304171')
    >>> url = 'https://github.com/mvdbeek/tools-devteam/raw/data_manager/data_managers/data_manager_fetch_genome_dbkeys_all_fasta/test-data/test.tar.gz'
    >>> params = {'param_dict': {'reference_source': {'user_url': url}}}
    >>> fh = download_from_url(params=params, tmp_dir=tmp_dir)[0][0]
    >>> assert fh.readline().startswith('>FBtr0304171')
    >>> url = 'https://github.com/mvdbeek/tools-devteam/raw/data_manager/data_managers/data_manager_fetch_genome_dbkeys_all_fasta/test-data/test.zip'
    >>> params = {'param_dict': {'reference_source': {'user_url': url}}}
    >>> fh = download_from_url(params=params, tmp_dir=tmp_dir)[0][0]
    >>> assert fh.readline().startswith('>FBtr0304171')
    >>> url = 'https://raw.githubusercontent.com/galaxyproject/tools-devteam/master/data_managers/data_manager_fetch_genome_dbkeys_all_fasta/test-data/phiX174.fasta'
    >>> params = {'param_dict': {'reference_source': {'user_url': url}}}
    >>> fh = download_from_url(params=params, tmp_dir=tmp_dir)[0][0]
    >>> assert fh.readline().startswith('>phiX174')
    """
    urls = filter(bool, [x.strip() for x in params['param_dict']['reference_source']['user_url'].split('\n')])
    return [get_stream_reader(urlopen(url), tmp_dir) for url in urls]


def download_from_history(params, tmp_dir, **kwds):
    # TODO: allow multiple FASTA input files
    input_filename = params['param_dict']['reference_source']['input_fasta']
    if isinstance(input_filename, list):
        fasta_readers = [get_stream_reader(open(filename, 'rb'), tmp_dir) for filename in input_filename]
    else:
        fasta_readers = get_stream_reader(open(input_filename), tmp_dir)
    return fasta_readers


def copy_from_directory(data_manager_dict, params, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, tmp_dir):
    input_filename = params['param_dict']['reference_source']['fasta_filename']
    create_symlink = params['param_dict']['reference_source']['create_symlink'] == 'create_symlink'
    if create_symlink:
        data_table_entries = _create_symlink(input_filename, target_directory, dbkey, dbkey_name, sequence_id, sequence_name)
    else:
        if isinstance(input_filename, list):
            fasta_readers = [get_stream_reader(open(filename, 'rb'), tmp_dir) for filename in input_filename]
        else:
            fasta_readers = get_stream_reader(open(input_filename), tmp_dir)
        return fasta_readers
    for data_table_name, data_table_entry in data_table_entries:
        if data_table_entry:
            _add_data_table_entry(data_manager_dict, data_table_entry, data_table_name)


def _add_data_table_entry(data_manager_dict, data_table_entry, data_table_name):
    data_manager_dict['data_tables'] = data_manager_dict.get('data_tables', {})
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get('all_fasta', [])
    data_manager_dict['data_tables'][data_table_name].append(data_table_entry)
    return data_manager_dict


def _stream_fasta_to_file(fasta_stream, target_directory, sequence_id, close_stream=True):
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join(target_directory, fasta_base_filename)
    with open(fasta_filename, 'wb+') as fasta_writer:

        if isinstance(fasta_stream, list) and len(fasta_stream) == 1:
            fasta_stream = fasta_stream[0]

        if isinstance(fasta_stream, list):
            last_char = None
            for fh in fasta_stream:
                if last_char not in [None, '\n', '\r', b'\n', b'\r']:
                    fasta_writer.write(b'\n')
                while True:
                    data = fh.read(CHUNK_SIZE)
                    if data:
                        fasta_writer.write(data)
                        last_char = data[-1]
                    else:
                        break
                if close_stream:
                    fh.close()
        else:
            while True:
                data = fasta_stream.read(CHUNK_SIZE)
                if data:
                    fasta_writer.write(data)
                else:
                    break
            if close_stream:
                fasta_stream.close()
    return fasta_filename


def compute_fasta_length(fasta_file, out_file, keep_first_word=False):
    infile = fasta_file
    out = open(out_file, 'w')

    fasta_title = ''
    seq_len = 0

    first_entry = True

    for line in open(infile):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line[0] == '>':
            if not first_entry:
                if keep_first_word:
                    fasta_title = fasta_title.split()[0]
                out.write("%s\t%d\n" % (fasta_title[1:], seq_len))
            else:
                first_entry = False
            fasta_title = line
            seq_len = 0
        else:
            seq_len += len(line)

    # last fasta-entry
    if keep_first_word:
        fasta_title = fasta_title.split()[0]
    out.write("%s\t%d\n" % (fasta_title[1:], seq_len))
    out.close()


def _create_symlink(input_filename, target_directory, dbkey, dbkey_name, sequence_id, sequence_name):
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join(target_directory, fasta_base_filename)
    os.symlink(input_filename, fasta_filename)

    dbkey_dict = None
    if dbkey_name:
        # do len calc here
        len_base_name = "%s.len" % (dbkey)
        compute_fasta_length(fasta_filename, os.path.join(target_directory, len_base_name), keep_first_word=True)
        dbkey_dict = dict(value=dbkey, name=dbkey_name, len_path=len_base_name)

    return [('__dbkeys__', dbkey_dict), ('all_fasta', dict(value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_filename))]


REFERENCE_SOURCE_TO_DOWNLOAD = dict(ucsc=download_from_ucsc, ncbi=download_from_ncbi, url=download_from_url, history=download_from_history,
                                    directory=copy_from_directory)
SORTING_METHODS = dict(as_is=_sort_fasta_as_is, lexicographical=_sort_fasta_lexicographical, gatk=_sort_fasta_gatk, custom=_sort_fasta_custom)


def main():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dbkey_description', dest='dbkey_description', action='store', type="string", default=None, help='dbkey_description')
    (options, args) = parser.parse_args()

    filename = args[0]

    params = loads(open(filename).read())
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory, mode=0o755)
    data_manager_dict = {}

    dbkey, dbkey_name, sequence_id, sequence_name = get_dbkey_dbname_id_name(params, dbkey_description=options.dbkey_description)

    if dbkey in [None, '', '?']:
        raise Exception('"%s" is not a valid dbkey. You must specify a valid dbkey.' % (dbkey))

    # Create a tmp_dir, in case a zip file needs to be uncompressed
    tmp_dir = tempfile.mkdtemp()
    # Fetch the FASTA
    try:
        reference_source = params['param_dict']['reference_source']['reference_source_selector']
        fasta_readers = REFERENCE_SOURCE_TO_DOWNLOAD[reference_source](data_manager_dict=data_manager_dict,
                                                                       params=params,
                                                                       target_directory=target_directory,
                                                                       dbkey=dbkey,
                                                                       dbkey_name=dbkey_name,
                                                                       sequence_id=sequence_id,
                                                                       sequence_name=sequence_name,
                                                                       tmp_dir=tmp_dir)
        if fasta_readers:
            add_fasta_to_table(data_manager_dict=data_manager_dict,
                               fasta_readers=fasta_readers,
                               target_directory=target_directory,
                               dbkey=dbkey,
                               dbkey_name=dbkey_name,
                               sequence_id=sequence_id,
                               sequence_name=sequence_name,
                               params=params)

    finally:
        cleanup_before_exit(tmp_dir)
    # save info to json file
    open(filename, 'wb').write(dumps(data_manager_dict).encode())


if __name__ == "__main__":
    main()
