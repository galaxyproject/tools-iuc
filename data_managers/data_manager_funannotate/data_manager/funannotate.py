#!/usr/bin/env python

import argparse
import json
import operator
import os
import subprocess
import sys
import tarfile
from datetime import datetime

import requests

# Some additional busco/orthodb10 datasets that can be added to funannotate db
# Will probably not be needed anymore in future versions of funannotate when it
# will use a recent busco version
BUSCO_10_DATASETS_URL = "https://busco-data.ezlab.org/v5/data/lineages/{dataset}"
BUSCO_10_DATASETS = [
    "acidobacteria_odb10.2020-03-06.tar.gz",
    "aconoidasida_odb10.2020-08-05.tar.gz",
    "actinobacteria_class_odb10.2021-02-23.tar.gz",
    "actinobacteria_phylum_odb10.2021-02-23.tar.gz",
    "actinopterygii_odb10.2021-02-19.tar.gz",
    "agaricales_odb10.2020-08-05.tar.gz",
    "agaricomycetes_odb10.2020-08-05.tar.gz",
    "alphabaculovirus_odb10.2020-11-26.tar.gz",
    "alphaherpesvirinae_odb10.2020-11-26.tar.gz",
    "alphaproteobacteria_odb10.2021-02-23.tar.gz",
    "alteromonadales_odb10.2021-02-23.tar.gz",
    "alveolata_odb10.2020-09-10.tar.gz",
    "apicomplexa_odb10.2020-09-10.tar.gz",
    "aquificae_odb10.2021-02-23.tar.gz",
    "arachnida_odb10.2020-08-05.tar.gz",
    "archaea_odb10.2021-02-23.tar.gz",
    "arthropoda_odb10.2020-09-10.tar.gz",
    "ascomycota_odb10.2020-09-10.tar.gz",
    "aves_odb10.2021-02-19.tar.gz",
    "aviadenovirus_odb10.2020-11-26.tar.gz",
    "bacillales_odb10.2021-02-23.tar.gz",
    "bacilli_odb10.2021-02-23.tar.gz",
    "bacteria_odb10.2020-03-06.tar.gz",
    "bacteroidales_odb10.2021-02-23.tar.gz",
    "bacteroidetes-chlorobi_group_odb10.2021-02-23.tar.gz",
    "bacteroidetes_odb10.2021-02-23.tar.gz",
    "bacteroidia_odb10.2021-02-23.tar.gz",
    "baculoviridae_odb10.2020-11-26.tar.gz",
    "basidiomycota_odb10.2020-09-10.tar.gz",
    "bclasvirinae_odb10.2020-11-26.tar.gz",
    "betabaculovirus_odb10.2020-11-26.tar.gz",
    "betaherpesvirinae_odb10.2020-11-26.tar.gz",
    "betaproteobacteria_odb10.2021-02-23.tar.gz",
    "boletales_odb10.2020-08-05.tar.gz",
    "brassicales_odb10.2020-08-05.tar.gz",
    "burkholderiales_odb10.2021-02-23.tar.gz",
    "campylobacterales_odb10.2020-03-06.tar.gz",
    "capnodiales_odb10.2020-08-05.tar.gz",
    "carnivora_odb10.2021-02-19.tar.gz",
    "cellvibrionales_odb10.2020-03-06.tar.gz",
    "cetartiodactyla_odb10.2021-02-19.tar.gz",
    "chaetothyriales_odb10.2020-08-05.tar.gz",
    "cheoctovirus_odb10.2020-11-26.tar.gz",
    "chlamydiae_odb10.2020-03-06.tar.gz",
    "chlorobi_odb10.2020-03-06.tar.gz",
    "chloroflexi_odb10.2020-03-06.tar.gz",
    "chlorophyta_odb10.2020-08-05.tar.gz",
    "chordopoxvirinae_odb10.2020-11-26.tar.gz",
    "chromatiales_odb10.2020-03-06.tar.gz",
    "chroococcales_odb10.2020-03-06.tar.gz",
    "clostridia_odb10.2020-03-06.tar.gz",
    "clostridiales_odb10.2020-03-06.tar.gz",
    "coccidia_odb10.2020-08-05.tar.gz",
    "coriobacteriales_odb10.2020-03-06.tar.gz",
    "coriobacteriia_odb10.2020-03-06.tar.gz",
    "corynebacteriales_odb10.2020-03-06.tar.gz",
    "cyanobacteria_odb10.2021-02-23.tar.gz",
    "cyprinodontiformes_odb10.2021-02-19.tar.gz",
    "cytophagales_odb10.2021-02-23.tar.gz",
    "cytophagia_odb10.2021-02-23.tar.gz",
    "delta-epsilon-subdivisions_odb10.2021-02-23.tar.gz",
    "deltaproteobacteria_odb10.2021-02-23.tar.gz",
    "desulfobacterales_odb10.2020-03-06.tar.gz",
    "desulfovibrionales_odb10.2021-02-23.tar.gz",
    "desulfurococcales_odb10.2021-02-23.tar.gz",
    "desulfuromonadales_odb10.2020-03-06.tar.gz",
    "diptera_odb10.2020-08-05.tar.gz",
    "dothideomycetes_odb10.2020-08-05.tar.gz",
    "embryophyta_odb10.2020-09-10.tar.gz",
    "endopterygota_odb10.2020-09-10.tar.gz",
    "enquatrovirus_odb10.2021-05-05.tar.gz",
    "enterobacterales_odb10.2021-02-23.tar.gz",
    "entomoplasmatales_odb10.2020-03-06.tar.gz",
    "epsilonproteobacteria_odb10.2020-03-06.tar.gz",
    "euarchontoglires_odb10.2021-02-19.tar.gz",
    "eudicots_odb10.2020-09-10.tar.gz",
    "euglenozoa_odb10.2020-08-05.tar.gz",
    "eukaryota_odb10.2020-09-10.tar.gz",
    "eurotiales_odb10.2020-08-05.tar.gz",
    "eurotiomycetes_odb10.2020-08-05.tar.gz",
    "euryarchaeota_odb10.2021-02-23.tar.gz",
    "eutheria_odb10.2021-02-19.tar.gz",
    "fabales_odb10.2020-08-05.tar.gz",
    "firmicutes_odb10.2021-02-23.tar.gz",
    "flavobacteriales_odb10.2021-02-23.tar.gz",
    "flavobacteriia_odb10.2021-02-23.tar.gz",
    "fromanvirus_odb10.2020-11-26.tar.gz",
    "fungi_odb10.2021-06-28.tar.gz",
    "fusobacteria_odb10.2020-03-06.tar.gz",
    "fusobacteriales_odb10.2020-03-06.tar.gz",
    "gammaherpesvirinae_odb10.2020-11-26.tar.gz",
    "gammaproteobacteria_odb10.2021-02-23.tar.gz",
    "glires_odb10.2021-02-19.tar.gz",
    "glomerellales_odb10.2020-08-05.tar.gz",
    "guernseyvirinae_odb10.2020-11-26.tar.gz",
    "halobacteria_odb10.2021-02-23.tar.gz",
    "halobacteriales_odb10.2021-02-23.tar.gz",
    "haloferacales_odb10.2021-02-23.tar.gz",
    "helotiales_odb10.2020-08-05.tar.gz",
    "hemiptera_odb10.2020-08-05.tar.gz",
    "herpesviridae_odb10.2020-11-26.tar.gz",
    "hymenoptera_odb10.2020-08-05.tar.gz",
    "hypocreales_odb10.2020-08-05.tar.gz",
    "insecta_odb10.2020-09-10.tar.gz",
    "iridoviridae_odb10.2020-11-26.tar.gz",
    "lactobacillales_odb10.2020-03-06.tar.gz",
    "laurasiatheria_odb10.2021-02-19.tar.gz",
    "legionellales_odb10.2020-03-06.tar.gz",
    "leotiomycetes_odb10.2020-08-05.tar.gz",
    "lepidoptera_odb10.2020-08-05.tar.gz",
    "liliopsida_odb10.2020-09-10.tar.gz",
    "mammalia_odb10.2021-02-19.tar.gz",
    "metazoa_odb10.2021-02-24.tar.gz",
    "methanobacteria_odb10.2021-02-23.tar.gz",
    "methanococcales_odb10.2021-02-23.tar.gz",
    "methanomicrobia_odb10.2021-02-23.tar.gz",
    "methanomicrobiales_odb10.2021-02-23.tar.gz",
    "micrococcales_odb10.2021-02-23.tar.gz",
    "microsporidia_odb10.2020-08-05.tar.gz",
    "mollicutes_odb10.2020-03-06.tar.gz",
    "mollusca_odb10.2020-08-05.tar.gz",
    "mucorales_odb10.2020-08-05.tar.gz",
    "mucoromycota_odb10.2020-08-05.tar.gz",
    "mycoplasmatales_odb10.2020-03-06.tar.gz",
    "natrialbales_odb10.2021-02-23.tar.gz",
    "neisseriales_odb10.2021-02-23.tar.gz",
    "nematoda_odb10.2020-08-05.tar.gz",
    "nitrosomonadales_odb10.2020-03-06.tar.gz",
    "nostocales_odb10.2020-03-06.tar.gz",
    "oceanospirillales_odb10.2020-03-06.tar.gz",
    "onygenales_odb10.2020-08-05.tar.gz",
    "oscillatoriales_odb10.2021-02-23.tar.gz",
    "pahexavirus_odb10.2020-11-26.tar.gz",
    "passeriformes_odb10.2021-02-19.tar.gz",
    "pasteurellales_odb10.2021-02-23.tar.gz",
    "peduovirus_odb10.2021-02-23.tar.gz",
    "planctomycetes_odb10.2020-03-06.tar.gz",
    "plasmodium_odb10.2020-08-05.tar.gz",
    "pleosporales_odb10.2020-08-05.tar.gz",
    "poales_odb10.2020-08-05.tar.gz",
    "polyporales_odb10.2020-08-05.tar.gz",
    "poxviridae_odb10.2020-11-26.tar.gz",
    "primates_odb10.2021-02-19.tar.gz",
    "propionibacteriales_odb10.2020-03-06.tar.gz",
    "proteobacteria_odb10.2021-02-23.tar.gz",
    "pseudomonadales_odb10.2020-03-06.tar.gz",
    "rhizobiales_odb10.2020-03-06.tar.gz",
    "rhizobium-agrobacterium_group_odb10.2020-03-06.tar.gz",
    "rhodobacterales_odb10.2021-02-23.tar.gz",
    "rhodospirillales_odb10.2020-03-06.tar.gz",
    "rickettsiales_odb10.2020-03-06.tar.gz",
    "rudiviridae_odb10.2020-11-26.tar.gz",
    "saccharomycetes_odb10.2020-08-05.tar.gz",
    "sauropsida_odb10.2021-02-19.tar.gz",
    "selenomonadales_odb10.2020-03-06.tar.gz",
    "simplexvirus_odb10.2020-11-26.tar.gz",
    "skunavirus_odb10.2020-11-26.tar.gz",
    "solanales_odb10.2020-08-05.tar.gz",
    "sordariomycetes_odb10.2020-08-05.tar.gz",
    "sphingobacteriia_odb10.2020-03-06.tar.gz",
    "sphingomonadales_odb10.2021-02-23.tar.gz",
    "spirochaetales_odb10.2020-03-06.tar.gz",
    "spirochaetes_odb10.2021-02-23.tar.gz",
    "spirochaetia_odb10.2021-02-23.tar.gz",
    "spounavirinae_odb10.2020-11-26.tar.gz",
    "stramenopiles_odb10.2020-08-05.tar.gz",
    "streptomycetales_odb10.2020-03-06.tar.gz",
    "streptosporangiales_odb10.2020-03-06.tar.gz",
    "sulfolobales_odb10.2021-02-23.tar.gz",
    "synechococcales_odb10.2020-03-06.tar.gz",
    "synergistetes_odb10.2020-03-06.tar.gz",
    "tenericutes_odb10.2020-03-06.tar.gz",
    "tequatrovirus_odb10.2020-11-26.tar.gz",
    "teseptimavirus_odb10.2020-11-26.tar.gz",
    "tetrapoda_odb10.2021-02-19.tar.gz",
    "tevenvirinae_odb10.2021-02-23.tar.gz",
    "thaumarchaeota_odb10.2021-02-23.tar.gz",
    "thermoanaerobacterales_odb10.2020-03-06.tar.gz",
    "thermoplasmata_odb10.2021-02-23.tar.gz",
    "thermoproteales_odb10.2021-02-23.tar.gz",
    "thermoprotei_odb10.2021-02-23.tar.gz",
    "thermotogae_odb10.2020-03-06.tar.gz",
    "thiotrichales_odb10.2020-03-06.tar.gz",
    "tissierellales_odb10.2020-03-06.tar.gz",
    "tissierellia_odb10.2020-03-06.tar.gz",
    "tremellomycetes_odb10.2020-08-05.tar.gz",
    "tunavirinae_odb10.2020-11-26.tar.gz",
    "varicellovirus_odb10.2020-11-26.tar.gz",
    "verrucomicrobia_odb10.2020-03-06.tar.gz",
    "vertebrata_odb10.2021-02-19.tar.gz",
    "vibrionales_odb10.2020-03-06.tar.gz",
    "viridiplantae_odb10.2020-09-10.tar.gz",
    "xanthomonadales_odb10.2020-03-06.tar.gz",
]


def download_file(url, dest):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--partial', dest='partial', action='store_true', help='Only download a small subset of data (for testing)')
    parser.add_argument('--wget', dest='wget', action='store_true', help='Download using wget (instead of urllib)')
    parser.add_argument("version_id")
    parser.add_argument("datatable_name")
    parser.add_argument("galaxy_datamanager_filename")
    args = parser.parse_args()

    with open(args.galaxy_datamanager_filename) as fh:
        config = json.load(fh)

    output_directory = config.get("output_data", [{}])[0].get("extra_files_path", None)
    data_manager_dict = {}
    data_manager_dict["data_tables"] = config.get("data_tables", {})
    data_manager_dict["data_tables"][args.datatable_name] = data_manager_dict[
        "data_tables"
    ].get(args.datatable_name, [])

    os.mkdir(output_directory)
    cmd_args = ['funannotate', 'setup', '-d', output_directory, '-b', 'all']
    if args.partial:
        cmd_args += ['-i', 'merops', '-b', 'eukaryota']
    if args.wget:
        cmd_args += ['--wget']
    proc = subprocess.Popen(args=cmd_args, shell=False, cwd=output_directory)
    return_code = proc.wait()
    if return_code:
        print("Error downloading Funannotate database.", file=sys.stderr)
        sys.exit(return_code)

    # Download newer busco datasets from orthodb 10
    if args.partial:
        BUSCO_10_DATASETS = BUSCO_10_DATASETS[:1]

    for busco_dataset in BUSCO_10_DATASETS:
        print("Downloading additional busco orthodb10 dataset %s" % busco_dataset)
        dest_tar = os.path.join(output_directory, busco_dataset)
        download_file(BUSCO_10_DATASETS_URL.format(dataset=busco_dataset), dest_tar)
        print("Extracting %s" % busco_dataset)
        tar = tarfile.open(dest_tar, "r:gz")
        tar.extractall(output_directory)
        tar.close()
        os.remove(dest_tar)

    version_id = datetime.today().strftime('%Y-%m-%d-%H%M%S')

    version = '1.0'

    data_manager_dict["data_tables"][args.datatable_name].append(
        dict(
            value=version_id,
            description="Funannotate database %s" % version_id,
            format_version=version,
            path=output_directory,
        )
    )

    data_manager_dict["data_tables"][args.datatable_name].sort(
        key=operator.itemgetter("value"), reverse=True
    )
    with open(args.galaxy_datamanager_filename, "w") as fh:
        json.dump(data_manager_dict, fh, indent=2, sort_keys=True)
