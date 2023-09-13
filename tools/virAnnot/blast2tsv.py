#!/usr/bin/env python3


# Name: blast2tsv
# Author(s): Sebastien Theil, Marie Lefebvre - INRAE
# Aims: Convert blast xml output to tsv and add taxonomy


import argparse
import logging as log

from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIXML
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def main():
    options = _set_options()
    _set_log_level(options.verbosity)
    hits = _read_xml(options)
    _write_csv(options, hits)


def _guess_database(accession):
    """Guess the correct database for querying based off the format of the accession"""
    database_mappings_refseq = {'AC_': 'nuccore', 'NC_': 'nuccore', 'NG_': 'nuccore',
                                'NT_': 'nuccore', 'NW_': 'nuccore', 'NZ_': 'nuccore',
                                'AP_': 'protein', 'NP_': 'protein', 'YP_': 'protein',
                                'XP_': 'protein', 'WP_': 'protein'}
    return database_mappings_refseq[accession[0:3]]


def _read_xml(options):
    """
    Parse XML blast results file
    Keep only the first hit
    """
    log.info("Read XML file.")
    results = open(options.xml_file, 'r')
    records = NCBIXML.parse(results)
    xml_results = {}
    for blast_record in records:
        for aln in blast_record.alignments:
            hit_count = 1
            for hit in aln.hsps:
                hsp = {}
                if hit_count == 1:
                    first_hit_frame = hit.frame[1] if len(hit.frame) > 0 else 0  # strand
                    cumul_hit_identity = hit.identities if hit.identities else 0
                    cumul_hit_score = hit.bits  # hit score
                    cumul_hit_evalue = hit.expect  # evalue
                    cumul_hit_length = hit.align_length if hit.align_length is not None else 0
                    hit_count = hit_count + 1
                else:
                    # all HSPs in different strand than 1st HSPs will be discarded.
                    if (first_hit_frame > 0 and hit.frame[1] > 0) or (first_hit_frame < 0 and hit.frame[1] < 0):
                        cumul_hit_identity = cumul_hit_identity + hit.identities
                        cumul_hit_length = cumul_hit_length + hit.align_length
                        cumul_hit_evalue = cumul_hit_evalue + hit.expect
                        cumul_hit_score = cumul_hit_score + hit.bits
                        hit_count = hit_count + 1
            if hit_count == 1:
                final_hit_count = hit_count
            elif hit_count > 1:
                final_hit_count = hit_count - 1
            hsp["evalue"] = cumul_hit_evalue / final_hit_count  # The smaller the E-value, the better the match
            hsp["query_id"] = blast_record.query_id
            hsp["query_length"] = blast_record.query_length  # length of the query
            hsp["accession"] = aln.accession.replace("ref|", "")
            hsp["description"] = aln.hit_def
            hsp["hit_length"] = aln.length  # length of the hit
            hsp["hsp_length"] = hit.align_length  # length of the hsp alignment
            hsp["queryOverlap"] = _get_overlap_value(options.algo, hsp, 'hsp', hsp["query_length"])[0]
            if cumul_hit_length == 0:
                hsp["percentIdentity"] = round(cumul_hit_identity, 1)  # identity percentage
            else:
                hsp["percentIdentity"] = round(cumul_hit_identity / cumul_hit_length * 100, 1)  # identity percentage
            hsp["score"] = cumul_hit_score  # The higher the bit-score, the better the sequence similarity
            hsp["num_hsps"] = final_hit_count
            hsp["hit_cumul_length"] = cumul_hit_length
            hsp["hitOverlap"] = _get_overlap_value(options.algo, hsp, 'hit', hsp["query_length"])[1]
            db = _guess_database(hsp["accession"])
            try:
                handle = Entrez.esummary(db=db, id=hsp["accession"])
                taxid = str(int(Entrez.read(handle)[0]['TaxId']))
                handle.close()
                log.info("Taxid found for " + hsp["accession"])
                lineage = ncbi.get_lineage(taxid)
                names = ncbi.get_taxid_translator(lineage)
                ordered = [names[tid] for tid in lineage]
                taxonomy = ordered[1:]
                hsp["tax_id"] = taxid
                hsp["taxonomy"] = ';'.join(taxonomy)
                hsp["organism"] = taxonomy[-1]
            except RuntimeError:
                hsp["tax_id"] = ""
                hsp["taxonomy"] = ""
                hsp["organism"] = ""
                log.warning("RuntimeError - Taxid not found for " + hsp["accession"])
            if hsp["evalue"] <= options.max_evalue and hsp["queryOverlap"] >= options.min_qov and \
                    hsp["hitOverlap"] >= options.min_hov and hsp["score"] >= options.min_score:
                xml_results[hsp["query_id"]] = hsp
            else:
                xml_results[hsp["query_id"]] = [hsp["query_length"]]

    return xml_results


def _get_overlap_value(algo, hsp, type, qlength):
    """
    Set hsp or hit overlap values for hit and query
    Return array [query_overlap, hit_overlap]
    """
    if type == 'hsp':
        q_align_len = qlength
        h_align_len = hsp["hsp_length"]
    else:
        q_align_len = qlength
        h_align_len = hsp["hit_cumul_length"]

    if algo == 'BLASTX':
        if q_align_len:
            query_overlap = (q_align_len * 3 / q_align_len) * 100
        if hsp["hit_length"]:
            hit_overlap = (h_align_len / hsp["hit_length"]) * 100
    elif algo == 'TBLASTN':
        if q_align_len:
            query_overlap = (q_align_len / q_align_len) * 100
        if hsp["hit_length"]:
            hit_overlap = (h_align_len * 3 / hsp["hit_length"]) * 100
    elif algo == 'TBLASTX':
        if q_align_len:
            query_overlap = (q_align_len * 3 / hsp["hsp_length"]) * 100
        if hsp["hit_length"]:
            hit_overlap = (h_align_len * 3 / hsp["hit_length"]) * 100
    else:
        if q_align_len:
            query_overlap = (q_align_len / q_align_len) * 100
        if hsp["hit_length"]:
            hit_overlap = (h_align_len / hsp["hit_length"]) * 100
    if query_overlap is None:
        query_overlap = 0
    if query_overlap > 100:
        query_overlap = 100
    if 'hit_overlap' not in locals():
        hit_overlap = 0
    if hit_overlap > 100:
        hit_overlap = 100

    return [round(query_overlap, 0), round(hit_overlap, 0)]


def _write_csv(options, hits):
    """
    Write output
    """
    # get a list of contig without corresponding number of mapped reads
    if options.rn_file is not None:
        with open(options.rn_file) as rn:
            rows = (line.split('\t') for line in rn)
            rn_list = {row[0]: row[1:] for row in rows}
    fasta = SeqIO.to_dict(SeqIO.parse(open(options.fasta_file), 'fasta'))
    headers = "#algo\tquery_id\tnb_reads\tquery_length\taccession\tdescription\torganism\tpercentIdentity\tnb_hsps\tqueryOverlap\thitOverlap\tevalue\tscore\ttax_id\ttaxonomy\tsequence\n"
    log.info("Write output file.")
    f = open(options.output, "w+")
    f.write(headers)
    for h in hits:
        if options.rn_file is not None:
            read_nb = ''.join(rn_list[h]).replace("\n", "")
        else:
            read_nb = ''
        if len(hits[h]) > 1:
            f.write(options.algo + "\t" + h + "\t" + read_nb + "\t" + str(hits[h]["query_length"]) + "\t")
            f.write(hits[h]["accession"] + "\t" + hits[h]["description"] + "\t")
            f.write(hits[h]["organism"] + "\t" + str(hits[h]["percentIdentity"]) + "\t")
            f.write(str(hits[h]["num_hsps"]) + "\t" + str(hits[h]["queryOverlap"]) + "\t")
            f.write(str(hits[h]["hitOverlap"]) + "\t" + str(hits[h]["evalue"]) + "\t")
            f.write(str(hits[h]["score"]) + "\t" + str(hits[h]["tax_id"]) + "\t")
            f.write(hits[h]["taxonomy"] + "\t" + str(fasta[h].seq))
            f.write("\n")
        else:
            f.write(options.algo + "\t" + h + "\t" + read_nb + "\t" + str(hits[h])[1:-1] + "\t")
            f.write("\n")
    f.close()


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--xml', help='XML files with results of blast', action='store', required=True, dest='xml_file')
    parser.add_argument('-rn', '--read-count', help='Tab-delimited file associating seqID with read number.', action='store', dest='rn_file')
    parser.add_argument('-c', '--contigs', help='FASTA file with contigs sequence.', action='store', required=True, dest='fasta_file')
    parser.add_argument('-me', '--max_evalue', help='Max evalue', action='store', type=float, default=0.0001, dest='max_evalue')
    parser.add_argument('-qov', '--min_query_overlap', help='Minimum query overlap', action='store', type=int, default=5, dest='min_qov')
    parser.add_argument('-mhov', '--min_hit_overlap', help='Minimum hit overlap', action='store', type=int, default=5, dest='min_hov')
    parser.add_argument('-s', '--min_score', help='Minimum score', action='store', type=int, default=30, dest='min_score')
    parser.add_argument('-a', '--algo', help='Blast type detection (BLASTN|BLASTP|BLASTX|TBLASTX|TBLASTN|DIAMONDX).', action='store', type=str, default='BLASTX', dest='algo')
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
