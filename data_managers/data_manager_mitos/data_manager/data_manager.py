import argparse
import json
import os
import shutil
import tarfile
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen

ZENODO = {
    "mitos": "2683856",
    "mitos2": "4284483"
}
NAMES = {
    "mitos1-refdata": "RefSeq39 + MiTFi tRNA models",
    "refseq39": "RefSeq39 (equivalent to MITOS1 data)",
    "refseq63m": "RefSeq63 Metazoa",
    "refseq63f": "RefSeq63 Fungi",
    "refseq63o": "RefSeq63 Opisthokonta",
    "refseq89m": "RefSeq89 Metazoa",
    "refseq89f": "RefSeq89 Fungi",
    "refseq89o": "RefSeq89 Opisthokonta"
}


def url_download(tpe, db, workdir):
    """
    download http://ab.inf.uni-tuebingen.de/data/software/megan6/download/FNAME
    to workdir
    and unzip

    return the name of the resulting dir
    """
    tarfname = os.path.join(workdir, db + ".tar.bz")
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    src = None
    dst = None
    try:
        req = Request("https://zenodo.org/record/{zenodoid}/files/{db}.tar.bz2?download=1".format(zenodoid=ZENODO[tpe], db=db))
        src = urlopen(req)
        with open(tarfname, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    finally:
        if src:
            src.close()
    with tarfile.open(tarfname, "r:bz2") as tar:
        dirname = tar.getnames()[0]
        tar.extractall(workdir)
    os.remove(tarfname)
    return dirname


def main(tpe, db, outjson):
    workdir = os.getcwd()

    path = url_download(tpe, db, workdir)

    data_manager_entry = {}
    data_manager_entry['value'] = db
    data_manager_entry['name'] = NAMES[db]
    data_manager_entry['type'] = tpe
    data_manager_entry['path'] = path
    data_manager_json = dict(data_tables=dict(mitos=data_manager_entry))

    with open(outjson) as fh:
        params = json.load(fh)
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)
    # output_path = os.path.abspath(os.path.join(os.getcwd(), 'mitos'))
    shutil.move(os.path.join(workdir, path), target_directory)
    with open(outjson, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--out', action='store', help='JSON filename')
    parser.add_argument('--type', action='store', help='mitos version')
    parser.add_argument('--db', action='store', help='db name')
    args = parser.parse_args()

    main(args.type, args.db, args.out)
