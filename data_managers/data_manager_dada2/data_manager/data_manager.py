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
    "silva_138.2_toGenus": "Silva version 138.2 to Genus (for short read 16S)",
    "silva_138.2_toSpecies": "Silva version 138.2 to Species (for long read 16S)",
    "silva_138": "Silva version 138",
    "silva_132": "Silva version 132",
    "silva_128": "Silva version 128",
    "rdp_19_toGenus": "RDP trainset 19 to Genus (for short read 16S)",
    "rdp_19_toSpecies": "RDP trainset 19 to Species (for long read 16S)",
    "rdp_16": "RDP trainset 16",
    "rdp_14": "RDP trainset 14",
    "greengenes_13.84": "GreenGenes version 13.84",
    "greengenes2_2024.09_toGenus": "GreenGenes2 release 2024.09 to Genus (for short read 16S)",
    "greengenes2_2024.09_toSpecies": "GreenGenes2 release 2024.09 to Species (for long read 16S)",
    "unite_8.0_fungi": "UNITE: General Fasta release 8.0 for Fungi",
    "unite_8.0_fungi_singletons": "UNITE: General Fasta release 8.0 for Fungi including global and 97% singletons",
    "RefSeq_RDP_2018_05": "NCBI RefSeq 16S rRNA database supplemented by RDP (05/2018)",
    "gtdb_2018_11": "GTDB: Genome Taxonomy Database (Bacteria &amp; Archaea) (11/2018)",
    "hitdb_1": "HitDB version 1 (Human InTestinal 16S rRNA)",
    "silva_euk_18S_132": "Silva version 132 Eukaryotic 18S",
    "PR2_4.11.1": "Protist Ribosomal Reference database (PR2) 4.11.1"
}

FILE2TAXURL = {
    # see https://github.com/benjjneb/dada2/issues/2053 for the distinction of the datasets
    "silva_138.2_toGenus": "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1",
    "silva_138.2_toSpecies": "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1",
    "silva_138": "https://zenodo.org/record/3731176/files/silva_nr_v138_train_set.fa.gz?download=1",
    "silva_132": "https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1",
    "silva_128": "https://zenodo.org/record/824551/files/silva_nr_v128_train_set.fa.gz?download=1",
    "rdp_19_toGenus": "https://zenodo.org/records/14168771/files/rdp_19_toGenus_trainset.fa.gz?download=1",
    "rdp_19_toSpecies": "https://zenodo.org/records/14168771/files/rdp_19_toSpecies_trainset.fa.gz?download=1",
    "rdp_16": "https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz?download=1",
    "rdp_14": "https://zenodo.org/record/158955/files/rdp_train_set_14.fa.gz?download=1",
    "unite_8.0_fungi": "https://s3.hpc.ut.ee/plutof-public/original/9f7b41c3-825b-4db8-9c52-74a4603a860a.zip",
    "unite_8.0_fungi_singletons": "https://s3.hpc.ut.ee/plutof-public/original/53dfc9ce-9cb5-4205-84bb-f47faff26462.zip",
    "greengenes_13.84": "https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1",
    "greengenes2_2024.09_toGenus": "https://zenodo.org/records/14169078/files/gg2_2024_09_toGenus_trainset.fa.gz?download=1",
    "greengenes2_2024.09_toSpecies": "https://zenodo.org/records/14169078/files/gg2_2024_09_toSpecies_trainset.fa.gz?download=1",
    "RefSeq_RDP_2018_05": "https://zenodo.org/record/2541239/files/RefSeq-RDP16S_v2_May2018.fa.gz?download=1",
    "gtdb_2018_11": "https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1",
    "hitdb_1": "https://zenodo.org/record/159205/files/hitdb_v1.00.fa.gz?download=1",
    "silva_euk_18S_132": "https://zenodo.org/record/1447330/files/silva_132.18s.99_rep_set.dada2.fa.gz?download=1",
    "PR2_4.11.1": "https://github.com/pr2database/pr2database/releases/download/4.11.1/pr2_version_4.11.1_dada2.fasta.gz"
}

FILE2SPECIESURL = {
    "silva_138.2": "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1",
    "silva_138": "https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1",
    "silva_132": "https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1",
    "silva_128": "https://zenodo.org/record/824551/files/silva_species_assignment_v128.fa.gz?download=1",
    # rdp_19 not available for assignSpecies https://github.com/benjjneb/dada2/issues/2053#issuecomment-2512185929
    "rdp_16": "https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz?download=1",
    "rdp_14": "https://zenodo.org/record/158955/files/rdp_species_assignment_14.fa.gz?download=1"
}

FILE2TAXLEVELS = {
    "PR2_4.11.1": "Kingdom,Supergroup,Division,Class,Order,Family,Genus,Species",
}


def url_download(url, fname, workdir):
    """
    download url to workdir/fname
    """
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

    # special treatment of UNITE DBs: they are zip files containing two fasta (xyz.fasta and developer/xyz.fasta)
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

    with open(outjson) as fh:
        params = json.load(fh)

    workdir = params['output_data'][0]['extra_files_path']
    os.mkdir(workdir)
    url_download(FILE2TAXURL[dataset], dataset + ".taxonomy", workdir)

    data_manager_json = {"data_tables": {}}
    data_manager_entry = {}
    data_manager_entry['value'] = dataset
    data_manager_entry['name'] = FILE2NAME[dataset]
    data_manager_entry['path'] = dataset + ".taxonomy"
    data_manager_entry['taxlevels'] = FILE2TAXLEVELS.get(dataset, DEFAULT_TAXLEVELS)
    data_manager_json["data_tables"]["dada2_taxonomy"] = data_manager_entry

    if FILE2SPECIESURL.get(dataset, False):
        url_download(FILE2SPECIESURL[dataset], dataset + ".species", workdir)
        data_manager_entry = {}
        data_manager_entry['value'] = dataset
        data_manager_entry['name'] = FILE2NAME[dataset]
        data_manager_entry['path'] = dataset + ".species"
        data_manager_json["data_tables"]["dada2_species"] = data_manager_entry

    with open(outjson, 'w') as fh:
        json.dump(data_manager_json, fh, sort_keys=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create data manager json.')
    parser.add_argument('--out', action='store', help='JSON filename')
    parser.add_argument('--dataset', action='store', help='Download data set name')
    args = parser.parse_args()

    remote_dataset(args.dataset, args.out)
