#!/usr/bin/env python3


# Name: virAnnot_otu
# Author: Marie Lefebvre - INRAE
# Reuirements: Ete3 toolkit and external apps
# Aims: Create viral OTUs based on RPS and Blast annotations


import argparse
import csv
import logging as log
import os
import random
import re

import pandas as pd
import xlsxwriter
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from ete3 import NodeStyle, SeqGroup, SeqMotifFace, Tree, TreeStyle


def main():
    """
    1 - retrieve info (sequence, query_id, taxo) from RPS file
    2 - align protein sequences of the same domain, calculate
    matrix of distances, generate trees
    3 - get statistics (read number) per otu
    4 - create HTML report
    """
    options = _set_options()
    _set_log_level(options.verbosity)
    hits_collection = _cut_sequence(options)
    _align_sequences(options, hits_collection)
    _get_stats(options, hits_collection)
    _create_html(options, hits_collection)


def _cut_sequence(options):
    """
    Retrieve viral hits and sequences from RPS files
    """
    log.info("Cut sequences")
    i = 0  # keep track of iterations over rps files to use the corresponding fasta file
    collection = {}
    options.rps.sort()
    for rps_file in options.rps:
        log.debug("Reading rps file " + str(rps_file))
        with open(rps_file[0], 'r') as rps_current_file:
            rps_reader = csv.reader(rps_current_file, delimiter='\t')
            headers = 0
            for row in rps_reader:
                if headers == 0:
                    # headers
                    headers += 1
                else:
                    if row[1] == "no_hit":
                        pass
                    else:
                        query_id = row[0]
                        cdd_id = row[2]
                        startQ = int(row[5])
                        endQ = int(row[6])
                        frame = float(row[7])
                        description = row[8]
                        superkingdom = row[9]
                        match = re.search("Viruses", superkingdom)
                        # if contig is viral then retrieve sequence
                        if match:
                            options.fasta.sort()
                            seq = _retrieve_fasta_seq(options.fasta[i][0], query_id)
                            seq_length = len(seq)
                            if endQ < seq_length:
                                seq = seq[startQ - 1:endQ]
                            else:
                                seq = seq[startQ - 1:seq_length]
                            if frame < 0:
                                seq = seq.reverse_complement()
                            prot = seq.translate()
                            if len(prot) >= options.min_protein_length:
                                log.debug("Add " + query_id + " to collection")
                                if cdd_id not in collection:
                                    collection[cdd_id] = {}
                                collection[cdd_id][query_id] = {}
                                collection[cdd_id][query_id]["nuccleotide"] = seq
                                collection[cdd_id][query_id]["protein"] = prot
                                collection[cdd_id][query_id]["full_description"] = description
                                if options.blast is not None:
                                    options.blast.sort()
                                    with open(options.blast[i][0], 'r') as blast_current_file:
                                        blast_reader = csv.reader(blast_current_file, delimiter='\t')
                                        for b_query in blast_reader:
                                            if b_query[1] == query_id:
                                                collection[cdd_id][query_id]["nb"] = b_query[2]
                                                if len(b_query) > 10:
                                                    collection[cdd_id][query_id]["taxonomy"] = b_query[14]
                                                else:
                                                    collection[cdd_id][query_id]["taxonomy"] = "Unknown"
                                            else:
                                                if "nb" not in collection[cdd_id][query_id]:
                                                    collection[cdd_id][query_id]["nb"] = 0
                                                if "taxonomy" not in collection[cdd_id][query_id]:
                                                    collection[cdd_id][query_id]["taxonomy"] = "Unknown"
                                else:
                                    log.info("No blast file")
                                    collection[cdd_id][query_id]["taxonomy"] = "Unknown"
                                    collection[cdd_id][query_id]["nb"] = 0

                                collection[cdd_id]["short_description"] = description.split(",")[0] + description.split(",")[1]  # keep pfamXXX and RdRp 1
                                collection[cdd_id]["full_description"] = description
        i += 1
    return collection


def _retrieve_fasta_seq(fasta_file, query_id):
    """
    From fasta file retrieve specific sequence with id
    """
    contigs_list = SeqIO.to_dict(SeqIO.parse(open(fasta_file), 'fasta'))
    seq = contigs_list[query_id].seq
    return seq


def _create_tree(tree, fasta, out, color):
    """
    Create phylogenic tree from multiple alignments
    """
    try:
        f = open(tree, 'r')
    except IOError:
        log.info("Unknown file: " + tree + ". You may have less than 2 sequences to align.")
        return

    line = ""
    for word in f:
        line += word.strip()

    f.close()
    seqs = SeqGroup(fasta, format="fasta")
    t = Tree(tree)
    ts = TreeStyle()
    ts.show_branch_length = True
    colors = _parse_color_file(color)
    node_names = t.get_leaf_names()
    for name in node_names:
        seq = seqs.get_seq(name)
        seqFace = SeqMotifFace(seq, seq_format="()")
        node = t.get_leaves_by_name(name)
        for i in range(0, len(node)):
            if name in colors:
                ns = NodeStyle()
                ns['bgcolor'] = colors[name]
                node[i].set_style(ns)
            node[i].add_face(seqFace, 0, 'aligned')

    t.render(out, tree_style=ts)


def _parse_color_file(file):
    fh = open(file)
    reader = csv.reader(fh, delimiter="\t")
    data = list(reader)
    colors = {}
    for i in range(0, len(data)):
        colors[data[i][0]] = data[i][1]

    return colors


def _align_sequences(options, hits_collection):
    """
    Align hit sequences with pfam reference
    """
    log.info("Align sequences")
    if not os.path.exists(options.output):
        os.mkdir(options.output)
    color_by_sample = {}
    for cdd_id in hits_collection:
        cdd_output = options.output + "/" + hits_collection[cdd_id]["short_description"].replace(" ", "_")
        if not os.path.exists(cdd_output):
            os.mkdir(cdd_output)
        if os.path.exists(cdd_output + "/seq_to_align.fasta"):
            os.remove(cdd_output + "/seq_to_align.fasta")
        file_seq_to_align = cdd_output + "/seq_to_align.fasta"
        file_color_config = cdd_output + "/color_config.txt"
        f = open(file_seq_to_align, "a")
        f_c = open(file_color_config, "w+")
        log.info("Writing to " + file_seq_to_align)
        count = 0  # count number of contig per domain
        for query_id in hits_collection[cdd_id]:
            if query_id not in ["short_description", "full_description"]:
                sample = query_id.split("_")[0]  # get sample from SAMPLE_IdCONTIG
                sample_color = "#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])
                # same color for each contig of the same sample
                if sample not in color_by_sample.keys():
                    color_by_sample[sample] = sample_color
                f.write(">" + query_id + "\n")
                f.write(str(hits_collection[cdd_id][query_id]["protein"]) + "\n")
                f_c.write(query_id + '\t' + color_by_sample[sample] + '\n')
                count += 1
        f.close()
        f_c.close()
        file_seq_aligned = cdd_output + '/seq_aligned.final_tree.fa'
        tree_file = cdd_output + '/tree.dnd'
        file_cluster = cdd_output + '/otu_cluster.csv'
        # create alignment for domain with more than 1 contigs
        if count > 1:
            log.info("Run clustal omega...")
            clustalo_cmd = ClustalOmegaCommandline("clustalo", infile=file_seq_to_align, outfile=file_seq_aligned,
                                                   guidetree_out=tree_file, seqtype="protein", force=True)
            log.debug(clustalo_cmd)
            stdout, stderr = clustalo_cmd()
            log.debug(stdout + stderr)

            # create tree plot with colors
            file_matrix = cdd_output + "/identity_matrix.csv"
            log.info("Create tree...")
            _create_tree(tree_file, file_seq_aligned, tree_file + '.png', file_color_config)
            _compute_pairwise_distance(options, file_seq_aligned, file_matrix, cdd_id)
            log.info("Retrieve OTUs...")
            # if os.path.exists(file_cluster):
            #     os.remove(file_cluster)
            otu_cmd = os.path.join(options.tool_path, 'seek_otu.R') + ' ' + file_matrix + ' ' + file_cluster + ' ' + str(options.perc)
            log.debug(otu_cmd)
            os.system(otu_cmd)
        # only one contig
        else:
            mv_cmd = 'cp ' + file_seq_to_align + ' ' + file_seq_aligned
            log.debug(mv_cmd)
            os.system(mv_cmd)

            f = open(file_cluster, "w+")
            f.write('OTU_1,1,' + list(hits_collection[cdd_id].keys())[0] + ',')
            f.close()


def _compute_pairwise_distance(options, file_seq_aligned, file_matrix, cdd_id):
    """
    Calculate paiwise distance between aligned protein sequences
    from a cdd_id
    """
    log.info("Compute pairwise distance of " + cdd_id)
    matrix = {}
    for k1 in SeqIO.parse(file_seq_aligned, "fasta"):
        row = []
        for k2 in SeqIO.parse(file_seq_aligned, "fasta"):
            identic = 0
            compared = 0
            keep_pos = 0
            for base in k1:
                base2 = k2[keep_pos]
                # mutation, next
                if base == 'X' or base2 == 'X':
                    keep_pos += 1
                    continue
                # gap in both sequences, next
                if base == '-' and base2 == '-':
                    keep_pos += 1
                    continue
                # gap in one of the sequence, next
                if base == '-' or base2 == '-':
                    keep_pos += 1
                    continue
                # identity
                if base == base2:
                    identic += 1
                compared += 1
                keep_pos += 1
            # set minimum overlap to 20
            if compared == 0 or compared < 20:
                percentIdentity = 0
            else:
                percentIdentity = (identic / compared) * 100
            row.append(percentIdentity)
        matrix[k1.id] = row
    log.debug("Write " + file_matrix)
    f = open(file_matrix, "w+")
    for row in matrix:
        f.write(row + ',' + ', '.join(map(str, matrix[row])) + "\n")
    f.close()


def _get_stats(options, hits_collection):
    """
    Retrieve annotation and number of read
    for    each OTUs
    """
    file_xlsx = options.output + '/otu_stats.xlsx'  # Create a workbook
    workbook = xlsxwriter.Workbook(file_xlsx)
    log.info("Writing stats to " + file_xlsx)
    for cdd_id in hits_collection:
        otu_collection = {}
        cdd_output = options.output + "/" + hits_collection[cdd_id]["short_description"].replace(" ", "_")
        worksheet = workbook.add_worksheet(hits_collection[cdd_id]["short_description"])  # add a worksheet
        file_cluster = cdd_output + '/otu_cluster.csv'
        with open(file_cluster, 'r') as clust:
            otu_reader = csv.reader(clust, delimiter=',')
            samples_list = []
            for row in otu_reader:
                contigs_list = row[2:len(row) - 1]  # remove last empty column
                otu_collection[row[0]] = {}  # key -> otu number
                otu_collection[row[0]]['contigs_list'] = contigs_list
                for contig in contigs_list:
                    sample = contig.split('_')[0]
                    samples_list.append(sample) if sample not in samples_list else samples_list
                    if sample not in otu_collection[row[0]]:
                        otu_collection[row[0]][sample] = {}
                        otu_collection[row[0]][sample][contig] = {}
                        # add read number of the contig and annotation
                        if 'nb' in hits_collection[cdd_id][contig]:
                            otu_collection[row[0]][sample][contig]['nb'] = hits_collection[cdd_id][contig]["nb"]
                        else:
                            otu_collection[row[0]][sample][contig]['nb'] = 0
                        if 'taxonomy' in hits_collection[cdd_id][contig]:
                            otu_collection[row[0]][sample][contig]['taxonomy'] = hits_collection[cdd_id][contig]["taxonomy"]
                        else:
                            otu_collection[row[0]][sample][contig]['taxonomy'] = 'unknown'
                    else:
                        otu_collection[row[0]][sample][contig] = {}
                        # add read number of the contig and annotation
                        if 'nb' in hits_collection[cdd_id][contig]:
                            otu_collection[row[0]][sample][contig]['nb'] = hits_collection[cdd_id][contig]["nb"]
                        else:
                            otu_collection[row[0]][sample][contig]['nb'] = 0
                        if 'taxonomy' in hits_collection[cdd_id][contig]:
                            otu_collection[row[0]][sample][contig]['taxonomy'] = hits_collection[cdd_id][contig]["taxonomy"]
                        else:
                            otu_collection[row[0]][sample][contig]['taxonomy'] = 'unknown'
                    if 'taxonomy' in hits_collection[cdd_id][contig]:
                        otu_collection[row[0]]['global_taxonomy'] = hits_collection[cdd_id][contig]["taxonomy"]
                    else:
                        otu_collection[row[0]]['global_taxonomy'] = 'unknown'

        # calculate total number of reads for each sample of each OTU
        for otu in otu_collection:
            for sample in otu_collection[otu]:
                if sample not in ['contigs_list', 'global_taxonomy']:
                    total_nb_read = 0
                    for contig in otu_collection[otu][sample]:
                        total_nb_read += int(otu_collection[otu][sample][contig]['nb'])
                    otu_collection[otu][sample]['total_nb_read'] = total_nb_read
        row = 0
        column = 0
        item = '#OTU_name'
        worksheet.write(row, column, item)
        for samp in samples_list:
            column += 1
            worksheet.write(row, column, samp)
        worksheet.write(row, column + 1, 'taxonomy')
        worksheet.write(row, column + 2, 'contigs_list')
        row = 1
        # column = 0
        for otu in otu_collection:
            if isinstance(otu_collection[otu], dict):
                column = 0
                worksheet.write(row, column, otu)
                # prepare table with 0 in each cells
                for sample in otu_collection[otu]:
                    column = 1
                    for samp in samples_list:
                        worksheet.write(row, column, 0)
                        column += 1
                # fill in table with nb of read for each sample and each OTU
                for sample in otu_collection[otu]:
                    column = 1
                    for samp in samples_list:
                        if samp == sample:
                            worksheet.write(row, column, otu_collection[otu][sample]['total_nb_read'])
                        column += 1
                worksheet.write(row, len(samples_list) + 1, otu_collection[otu]['global_taxonomy'].replace(';', ' '))
                worksheet.write(row, len(samples_list) + 2, ",".join(otu_collection[otu]['contigs_list']))
                row += 1
    workbook.close()
    read_file = pd.ExcelFile(file_xlsx)
    for sheet in read_file.sheet_names:
        cluster_nb_reads_file = options.output + "/" + sheet.replace(" ", "_") + "/cluster_nb_reads_files.tab"
        data_xls = pd.read_excel(file_xlsx, sheet, dtype=str, index_col=None)
        data_xls.to_csv(cluster_nb_reads_file, encoding='utf-8', index=False, sep='\t')


def _create_html(options, hits_collection):
    """
    Create HTML file with all results
    """
    # create mapping file with all informations to use to create HTML report
    map_file_path = options.output + "/map.txt"
    if os.path.exists(map_file_path):
        os.remove(map_file_path)

    map_file = open(map_file_path, "w+")
    headers = ['#cdd_id', 'align_files', 'tree_files', 'cluster_files', 'cluster_nb_reads_files', 'pairwise_files', 'description', 'full_description\n']
    map_file.write("\t".join(headers))
    for cdd_id in hits_collection:
        cdd_output = hits_collection[cdd_id]["short_description"].replace(" ", "_")
        short_description = cdd_output
        file_seq_aligned = cdd_output + '/seq_aligned.final_tree.fa'
        tree_file = cdd_output + '/tree.dnd.png'
        file_cluster = cdd_output + '/otu_cluster.csv'
        file_matrix = cdd_output + "/identity_matrix.csv"
        cluster_nb_reads_files = cdd_output + "/cluster_nb_reads_files.tab"
        map_file.write(cdd_id + "\t" + file_seq_aligned + "\t" + tree_file + "\t")
        map_file.write(file_cluster + "\t" + cluster_nb_reads_files + "\t" + file_matrix + "\t")
        map_file.write(short_description + "\t" + hits_collection[cdd_id]["full_description"] + "\n")
    map_file.close()
    log.info("Writing HTML report")
    html_cmd = os.path.join(options.tool_path, 'rps2tree_html.py') + ' -m ' + map_file_path + ' -o ' + options.output
    log.debug(html_cmd)
    os.system(html_cmd)


def _set_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--blast', help='TAB blast file from blast2ecsv module.', action='append', required=False, dest='blast', nargs='+')
    parser.add_argument('-r', '--rps', help='TAB rpsblast file from rps2ecsv module.', action='append', required=True, dest='rps', nargs='+')
    parser.add_argument('-f', '--fasta', help='FASTA file with contigs', action='append', required=True, dest='fasta', nargs='+')
    parser.add_argument('-p', '--percentage', help='Percentage similarity threshold for OTUs cutoff.', action='store', type=int, default=90, dest='perc')
    parser.add_argument('-vp', '--viral_portion', help='Minimun portion of viral sequences in RPS domain to be included.', action='store', type=float, default=0.3, dest='viral_portion')
    parser.add_argument('-mpl', '--min_protein_length', help='Minimum query protein length.', action='store', type=int, default=100, dest='min_protein_length')
    parser.add_argument('-tp', '--tool_path', help='Path to otu_seek.R', action='store', type=str, default='./', dest='tool_path')
    parser.add_argument('-o', '--out', help='The output directory', action='store', type=str, default='./Rps2tree_OTU', dest='output')
    parser.add_argument('-rgb', '--rgb-conf', help='Color palette for contigs coloration', action='store', type=str, default='rgb.txt', dest='file_rgb')
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
