import argparse
import json
import os
try:
    # For Python 3.0 and later
    from urllib.request import Request, urlopen
except ImportError:
    # Fall back to Python 2 imports
    from urllib2 import Request, urlopen

DEFAULT_TAXLEVELS = "Kingdom,Phylum,Class,Order,Family,Genus,Species"

FILE2NAME = {
    "silva_132": "Silva version 132",
    "silva_128": "Silva version 128",
    "rdp_16": "RDP trainset 16",
    "rdp_14": "RDP trainset 14",
    "greengenes_13.84": "GreenGenes version 13.84",
    "unite_8.0_fungi": "UNITE: General Fasta release 8.0 for Fungi",
    "unite_8.0_fungi_singletons": "UNITE: General Fasta release 8.0 for Fungi including global and 97% singletons",
    "RefSeq_RDP_2018_05": "NCBI RefSeq 16S rRNA database supplemented by RDP (05/2018)",
    "gtdb_2018_11": "GTDB: Genome Taxonomy Database (Bacteria &amp; Archaea) (11/2018)",
    "hitdb_1": "HitDB version 1 (Human InTestinal 16S rRNA)",
    "silva_euk_18S_132": "Silva version 132 Eukaryotic 18S",
    "PR2_4.11.1": "Protist Ribosomal Reference database (PR2) 4.11.1"
}

FILE2TAXURL = {
    "silva_132": "https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1",
    "silva_128": "https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz?download=1",
    "rdp_16": "https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz?download=1",
    "rdp_14": "https://zenodo.org/record/158955/files/rdp_train_set_14.fa.gz?download=1",
    "unite_8.0_fungi": "https://files.plutof.ut.ee/public/orig/EB/0C/EB0CCB3A871B77EA75E472D13926271076904A588D2E1C1EA5AFCF7397D48378.zip",
    "unite_8.0_fungi_singletons": "https://files.plutof.ut.ee/doi/06/A2/06A2C86256EED64085670EB0C54B7115F6DAC8F311C656A9CB33E386CFABA0D0.zip",
    "greengenes_13.84": "https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1",
    "RefSeq_RDP_2018_05": "https://zenodo.org/record/2541239/files/RefSeq-RDP16S_v2_May2018.fa.gz?download=1",
    "gtdb_2018_11": "https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1",
    "hitdb_1": "https://zenodo.org/record/159205/files/hitdb_v1.00.fa.gz?download=1",
    "silva_euk_18S_132": "https://zenodo.org/record/1447330/files/silva_132.18s.99_rep_set.dada2.fa.gz?download=1",
    "PR2_4.11.1": "https://github.com/pr2database/pr2database/releases/download/4.11.1/pr2_version_4.11.1_dada2.fasta.gz"
}

FILE2SPECIESURL = {
    "silva_132": "https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1",
    "silva_128": "https://zenodo.org/record/824551/files/silva_species_assignment_v128.fa.gz?download=1",
    "rdp_16": "https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz?download=1",
    "rdp_14": "https://zenodo.org/record/158955/files/rdp_species_assignment_14.fa.gz?download=1"
}

FILE2TAXLEVELS = {
    "PR2_4.11.1": "Kingdom,Supergroup,Division,Class,Order,Family,Genus,Species"
}


def url_download(url, fname, workdir):
    """
    download url to workdir/fname
    """
    import logging
    logging.error("DL %s"%url)
    file_path = os.path.join(workdir, fname)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    src = None
    dst = None
    try:
        req = Request(url)
        src = urlopen(req)
        with open(file_path, 'wb') as dst:
            while True:
                chunk = src.read(2**10)
                if chunk:
                    dst.write(chunk)
                else:
                    break
    finally:
        if src:
            src.close()

#   special treatment of UNITE DBs: they are zip files containing two fasta (xyz.fasta and developer/xyz.fasta)
    if fname.startswith("unite"):
        import glob
        import gzip
        import shutil
        import zipfile
        # unzip download
        zip_ref = zipfile.ZipFile(file_path, 'r')
        zip_ref.extractall(workdir)
        zip_ref.close()
        # gzip top level fasta file
        fastas = glob.glob("%s/*fasta" % workdir)
        if len(fastas) != 1:
            msg = "UNITE download %s contained %d fasta file(s): %s" % (url, len(fastas), " ".join(fastas))
            raise Exception(msg)
        with open(fastas[0], 'rb') as f_in:
            with gzip.open(file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


def remote_dataset(dataset, outjson):

    with open(outjson) as jf:
        params = json.loads(jf.read())

    workdir = params['output_data'][0]['extra_files_path']
    os.mkdir(workdir)
    url_download( FILE2TAXURL[dataset], dataset + ".taxonomy", workdir)

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = dataset
    data_manager_entry['name'] = FILE2NAME[dataset]
    data_manager_entry['path'] = dataset + ".taxonomy"
    data_manager_entry['taxlevels'] = FILE2TAXLEVELS.get(dataset, DEFAULT_TAXLEVELS)
    data_manager_json["data_tables"]["dada2_taxonomy"] = data_manager_entry

    if FILE2SPECIESURL.get(dataset, False ):
        url_download( FILE2SPECIESURL[dataset], dataset + ".species", workdir)
        data_manager_entry = {}
        data_manager_entry['value'] = dataset
        data_manager_entry['name'] = FILE2NAME[dataset]
        data_manager_entry['path'] = dataset + ".species"
        data_manager_json["data_tables"]["dada2_species"] = data_manager_entry

    with open(outjson, 'w') as jf:
        jf.write(json.dumps(data_manager_json))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--out', action='store', help='JSON filename')
    parser.add_argument('--dataset', action='store', help='Download data set name')
    args = parser.parse_args()

    remote_dataset(args.dataset, args.out)
