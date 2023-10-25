#!/usr/bin/env python3


# Name: rps2ecsv
# Author: Marie Lefebvre - INRAE
# Aims: Convert rpsblast xml output to csv and add taxonomy


import argparse
import logging as log
import json
from urllib import request
from urllib.error import URLError, HTTPError

from Bio.Blast import NCBIXML
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def main():
    options = _set_options()
    _set_log_level(options.verbosity)
    hits = _read_xml(options)
    _write_tsv(options, hits)


def _read_xml(options):
    """
    Parse XML RPSblast results file
    """
    log.info("Read XML file " + options.xml_file)
    xml = open(options.xml_file, 'r')
    records = NCBIXML.parse(xml)
    xml_results = {}
    for blast_record in records:
        for aln in blast_record.alignments:
            for hit in aln.hsps:
                hsp = {}
                hit_evalue = hit.expect
                if hit_evalue > options.max_evalue:
                    continue
                hit_frame = hit.frame[0]  # frame
                hit_evalue = hit.expect  # evalue
                hit_startQ = hit.query_start
                hit_endQ = hit.query_end
            hsp["frame"] = hit_frame
            hsp["evalue"] = hit_evalue
            hsp["startQ"] = hit_startQ
            hsp["endQ"] = hit_endQ
            hsp["query_id"] = blast_record.query_id
            hsp["cdd_id"] = aln.hit_def.split(",")[0]
            hsp["hit_id"] = aln.hit_id
            hsp["query_length"] = blast_record.query_length  # length of the query
            hsp["description"] = aln.hit_def
            hsp["accession"] = aln.accession
            hsp["pfam_id"] = hsp["description"].split(",")[0].replace("pfam", "PF")
            log.info("Requeting Interpro for " + hsp["pfam_id"])
            url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/" + hsp["pfam_id"] + "/taxonomy/uniprot/"
            req = request.Request(url)
            try:
                response = request.urlopen(req)
            except HTTPError as e:
                log.debug('Http error for interpro: ', e.code)
            except URLError as e:
                log.debug('Url error for interpro: ', e.reason)
            else:
                encoded_response = response.read()
                decoded_response = encoded_response.decode()
                payload = json.loads(decoded_response)
                kingdoms = []
                for item in payload["taxonomy_subset"]:
                    lineage_string = item["lineage"]
                    lineage = [int(i) for i in lineage_string]
                    translation = ncbi.get_taxid_translator(lineage)
                    names = list(translation.values())
                    taxonomy = names[1:]  # remove 'root' at the begining
                    kingdoms.append(taxonomy[0])
                frequency = {kingdom: kingdoms.count(kingdom) for kingdom in kingdoms}  # {'Pseudomonadota': 9, 'cellular organisms': 4}
                sorted_freq = dict(sorted(frequency.items(), key=lambda x: x[1], reverse=True))
                concat_freq = ";".join("{}({})".format(k, v) for k, v in sorted_freq.items())
                hsp["taxonomy"] = concat_freq
            xml_results[hsp["query_id"]] = hsp
    return xml_results


def _write_tsv(options, hits):
    """
    Write output
    """
    log.info("Write output file " + options.output)
    headers = "#query_id\tquery_length\tcdd_id\thit_id\tevalue\tstartQ\tendQ\tframe\tdescription\tsuperkingdom\n"
    f = open(options.output, "w+")
    f.write(headers)
    for h in hits:
        f.write(h + "\t" + str(hits[h]["query_length"]) + "\t")
        f.write(hits[h]["cdd_id"] + "\t" + hits[h]["hit_id"] + "\t" + str(hits[h]["evalue"]) + "\t")
        f.write(str(hits[h]["startQ"]) + "\t" + str(hits[h]["endQ"]) + "\t" + str(hits[h]["frame"]) + "\t")
        f.write(hits[h]["description"] + "\t" + hits[h]["taxonomy"])
        f.write("\n")
    f.close()


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xml', help='XML files with results of blast', action='store', required=True, dest='xml_file')
    parser.add_argument('-e', '--max_evalue', help='Max evalue', action='store', type=float, default=0.0001, dest='max_evalue')
    parser.add_argument('-o', '--out', help='The output file (.tab).', action='store', type=str, default='./rps2tsv_output.tab', dest='output')
    parser.add_argument('-v', '--verbosity', help='Verbose level', action='store', type=int, choices=[1, 2, 3, 4], default=1)
    args = parser.parse_args()
    return args


def _set_log_level(verbosity):
    if verbosity == 1:
        log_format = '%(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.INFO, format=log_format)
    elif verbosity == 3:
        log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.DEBUG, format=log_format)


if __name__ == "__main__":
    main()
