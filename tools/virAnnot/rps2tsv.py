#!/usr/bin/env python3


# Name: rps2ecsv
# Author: Marie Lefebvre - INRAE
# Aims: Convert rpsblast xml output to csv and add taxonomy


import argparse
import gzip
import os
import logging as log
import subprocess
from Bio.Blast import NCBIXML


def main():
    options = _set_options()
    _set_log_level(options.verbosity)
    _get_db()
    hits = _read_xml(options)
    _write_tsv(options, hits)


def _get_db():
    try:
        from urllib import urlretrieve
    except ImportError:
        from urllib.request import urlretrieve
    if not os.path.exists("pfamA_tax_depth.txt"):
        log.info('Downloading pfamA_tax_depth.txt.gz from NCBI FTP site (via HTTP)...')
        urlretrieve("http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA_tax_depth.txt.gz", 'pfamA_tax_depth.txt.gz')
        log.info('Done. Parsing...')
        input = gzip.GzipFile("pfamA_tax_depth.txt.gz", 'rb')
        gz_content = input.read()
        input.close()
        output = open("pfamA_tax_depth.txt", "wb")
        output.write(gz_content)
        output.close()


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
            hsp["pfam_id"] = hsp["description"].split(",")[0]
            cmd = "grep " + hsp["pfam_id"].replace("pfam", "PF") + " " + options.pfamA_tax_depth
            log.debug(cmd)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            out = str(out).replace('b\'', '').replace('\'', '').split("\\n")
            hsp["taxonomy"] = ""
            for pf in out[:-1]:  # loop except last item which is empty
                pf = pf.split("\\t")
                if int(pf[2]) > 9:  # count
                    hsp["taxonomy"] += pf[1] + "(" + pf[2] + ");"
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
    parser.add_argument('-pf', '--pfam', help='pfamA_tax_depth.txt file http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/database_files/',
                        action='store', default='pfamA_tax_depth.txt', dest='pfamA_tax_depth')
    parser.add_argument('-o', '--out', help='The output file (.csv).', action='store', type=str, default='./', dest='output')
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
