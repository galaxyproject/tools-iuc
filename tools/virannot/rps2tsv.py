#!/usr/bin/env python3


# Name: rps2ecsv
# Author: Marie Lefebvre - INRAE
# Aims: Convert rpsblast xml output to csv and add taxonomy

"""Module which converts rpsblast xml output to tsv and add taxonomy"""

import argparse
import json
import logging as log
from urllib import request
from urllib.error import HTTPError, URLError

from Bio.Blast import NCBIXML
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def main():
    """
    Main function
    """
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
                hit_identity = hit.identities
                hit_aln_length = hit.align_length
                pident = "%0.3f" % (100 * float(hit_identity) / float(hit_aln_length))
            if float(pident) < 0.1:
                continue
            hsp["pident"] = pident
            hsp["frame"] = hit_frame
            hsp["evalue"] = hit_evalue
            hsp["startQ"] = hit_startQ
            hsp["endQ"] = hit_endQ
            hsp["query_id"] = blast_record.query
            hsp["cdd_id"] = aln.hit_def.split(",")[0]
            hsp["hit_id"] = aln.hit_id
            hsp["query_length"] = blast_record.query_length  # length of the query
            hsp["description"] = aln.hit_def
            hsp["accession"] = aln.accession
            hsp["pfam_id"] = hsp["description"].split(",")[0].replace("pfam", "PF")
            log.info("Requeting Interpro for " + hsp["pfam_id"])
            url = "https://www.ebi.ac.uk/interpro/api/taxonomy/uniprot/entry/pfam/" + hsp["pfam_id"]
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
                for item in payload["results"][:6]:
                    if item["metadata"]["parent"] is not None:
                        lineage_parent = item["metadata"]["parent"]
                        translation = ncbi.get_taxid_translator([int(lineage_parent)])
                        names = list(translation.values())
                        if len(names) > 0:
                            if names[0] == "root":
                                taxonomy = names[1:]  # remove 'root' at the begining
                            else:
                                taxonomy = names
                        else:
                            taxonomy = names
                        if len(taxonomy) != 0:
                            kingdoms.append(taxonomy[0])
                # {'Pseudomonadota': 9, 'cellular organisms': 4}
                frequency = {kingdom: kingdoms.count(kingdom) for kingdom in kingdoms}
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
    headers = "#query_id\tquery_length\tcdd_id\thit_id\tevalue\tstartQ\tendQ\tframe\tdescription\tsuperkingdom\tpident\n"
    f = open(options.output, "w+")
    f.write(headers)
    for h in hits:
        f.write(h + "\t" + str(hits[h]["query_length"]) + "\t")
        f.write(hits[h]["cdd_id"] + "\t" + hits[h]["hit_id"] + "\t" + str(hits[h]["evalue"]) + "\t")
        f.write(str(hits[h]["startQ"]) + "\t" + str(hits[h]["endQ"]) + "\t"
                + str(hits[h]["frame"]) + "\t")
        f.write(hits[h]["description"] + "\t" + hits[h]["taxonomy"] + "\t" + hits[h]["pident"])
        f.write("\n")
    f.close()


def _set_options():
    """
    Script parameters
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xml', help='XML files with results of blast', action='store',
                        required=True, dest='xml_file')
    parser.add_argument('-e', '--max_evalue', help='Max evalue', action='store',
                        type=float, default=0.0001, dest='max_evalue')
    parser.add_argument('-o', '--out', help='The output file (.tab).', action='store',
                        type=str, default='./rps2tsv_output.tab', dest='output')
    parser.add_argument('-v', '--verbosity', help='Verbose level', action='store',
                        type=int, choices=[1, 2, 3, 4], default=1)
    args = parser.parse_args()
    return args


def _set_log_level(verbosity):
    """
    Debbug
    """
    if verbosity == 1:
        log_format = '%(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.INFO, format=log_format)
    elif verbosity == 3:
        log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
        log.basicConfig(level=log.DEBUG, format=log_format)


if __name__ == "__main__":
    main()
