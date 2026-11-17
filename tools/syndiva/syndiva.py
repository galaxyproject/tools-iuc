#!/usr/bin/env python
# title : syndiva.py
# description : This script will analyze fasta files, look for restriction sites,
# cut the sequences around the restriction sites,
# translate the nucleic sequences into amino acids sequences.
# author : Fabienne Wong Jun Tai and Benjamin Dartigues
# creation date : 20121107
# version : 1.0 - revised November 2012
# version : 1.1 - revised March 2022
# usage : python syndiva.py -i file.fasta -o /output/dir/ -p pattern -5 seq_restric_5'-3 seq_restric_3'
# notes :
# # python_version :3.7.11
# # biopython_max_version  :1.72
# ==============================================================================
import math
import re
import subprocess
import sys

import matplotlib
import numpy
from args import Args
from args import get_os_path_join, get_os_path_name
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SubsMat import MatrixInfo

matplotlib.use('Agg')
from matplotlib import pyplot as plot  # noqa: I202,E402


args = Args()
# Variables initialization
directory = args.output_dir
mcl_file = get_os_path_join(directory, "mcl.in")
mcl_output = get_os_path_join(directory, "mcl.out")
html_file = get_os_path_join(directory, "syndiva_report.html")
graph_pic = get_os_path_join(directory, "distri.png")
input_file = get_os_path_name(args.input)
site_res_5 = args.site_res_5
site_res_3 = args.site_res_3
tag = {'mut': [], 'ok_stop_ext': [], 'stop': [], 'no_restric': [], 'no_multiple': [], 'amber': []}
all_seq = []
all_seq_fasta = {}  # dictionnary that will store information about all the sequences
good_seq = {}  # dictionnary that will store information about the valid sequences
identical_clones = {}
var_seq_common = {}  # dictionnary that will store the number of sequences that share the same variable parts
align_scores = []
nb_var_part = 0


def get_identity(str1, str2):
    if len(str2) > len(str1):
        return (len(str2) - len([i for i in range(len(str1)) if str1[i] != str2[i]])) / len(str2)
    else:
        return (len(str1) - len([i for i in range(len(str1)) if str1[i] != str2[i]])) / len(str1)


def reverse_complement(_seq):
    return str(Seq(_seq).reverse_complement())


def generate_aln(seq_dic, ids):  # sourcery skip: use-join
    # Multiple Sequence Alignment via ClustalO
    _input = ''
    for sequence_id in ids:
        _input += '>%s\n%s\n' % (sequence_id, re.sub("(.{80})", "\\1\n", seq_dic[sequence_id]['prot'], re.DOTALL))
    p = subprocess.Popen(["clustalo", "-i", "-", "--outfmt", "clu"], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
    aln_out, aln_err = p.communicate(input=_input)
    return aln_out


def report_html(_html_file, _tag, _all_seq, _good_seq, _all_seq_fasta, _identical_clones, _nb_var_part, _var_seq_common, _align_scores, _args):
    # Generate the html file for the report
    _all_seq.sort()
    for key in _tag.keys():
        _tag[key].sort()
    _good_seq = dict(sorted(_good_seq.items()))
    good_ids = _good_seq.keys()
    w = open(_html_file, 'w')
    w.write(
        '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN""http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"><html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" '
        'lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=utf-8" /><title>SynDivA Report</title><link '
        'href="http://twitter.github.com/bootstrap/assets/css/bootstrap.css" rel="stylesheet" /><style type="text/css">body {padding-top: 40px;}.subhead {padding: 40px '
        '0;}.subhead h1 {font-size: 60px;}.fasta {   font-family: Monaco, Menlo, Consolas, "Courier New", monospace;   font-size: 12px;}code.grey{color: '
        '#636D71;}</style></head><body><a id="top"></a><div class="navbar navbar-fixed-top"><div class="navbar-inner"><div class="container"><a class="brand" href="#top">SynDivA '
        'Report</a><div class="nav-collapse collapse"><ul class="nav"><li><a href="#input">Input data</a></li><li><a href="#analysis">Sequences analysis</a></li><li><a '
        'href="#variable">Variable regions analysis</a></li><li><a href="#cluster">Clustering</a></li><li><a href="#stat">Statistics</a></li><li><a '
        'href="#annex">Annex</a></li></ul></div></div></div></div><div class="container-fluid"><header class="subhead"><h1>SynDivA Report</h1></header><div '
        'class="page-header"><a id="input"></a><h2>Input data</h2></div>')

    # Input data
    w.write(
        '<p>Input file:<br/><code class="grey">%s</code></p><p>Number of sequences in input file:<br/><code class="grey">%d</code></p><p>Pattern of the sequence bank:<br/><code '
        'class="grey">%s</code></p><p>5\' restriction site:<br/><code class="grey">%s</code></p><p>3\' restriction site:<br/><code class="grey">%s</code></p>' % (
            input_file, len(_all_seq), _args.pattern, _args.site_res_5, _args.site_res_3))

    # Sequence analysis
    w.write(
        '<div class="page-header"><a id="analysis"></a><h2>Sequences analysis</h2></div><p>Caption:</p><ul><li class="text-success">Valid sequences that will be part of the next '
        'analysis </li><li class="text-warning">Good sequences but will not be part of the next analysis</li><li class="text-error">Rejected sequences</li></ul><table '
        'class="table table-striped table-bordered"><tr><th class="text-error">Absence of restriction sites</th><th class="text-error">Incorrect number of nucleotides between '
        'the restriction sites</th><th class="text-error">Stop codon <u>inside</u> the area of interest</th><th class="text-warning">Mutation in the conserved regions</th><th '
        'class="text-success">Valid sequences</th><th>Amber codon in the sequence (<u>inside</u> the area of interest)</th></tr>')
    w.write(
        '<tr><td class="text-error">%d sequence(s) (%.2f%%)</td><td class="text-error">%d sequence(s) (%.2f%%)</td><td class="text-error">%d sequence(s) (%.2f%%)</td><td '
        'class="text-warning">%d sequence(s) (%.2f%%)</td><td class="text-success">%d sequence(s) (%.2f%%)</td><td>%d sequence(s)</td></tr>' % (
            len(_tag['no_restric']), float(len(_tag['no_restric'])) / float(len(_all_seq)) * 100, len(_tag['no_multiple']), float(len(_tag['no_multiple'])) / float(len(_all_seq)) * 100, len(_tag['stop']),
            float(len(_tag['stop'])) / float(len(_all_seq)) * 100, len(_tag['mut']), float(len(_tag['mut'])) / float(len(_all_seq)) * 100, len(good_ids),
            float(len(good_ids)) / float(len(_all_seq)) * 100,
            len(_tag['amber'])))
    w.write(
        '<tr><td class="text-error">%s</td><td class="text-error">%s</td><td class="text-error">%s</td><td class="text-warning">%s</td><td '
        'class="text-success">%s</td><td>%s</td></tr></table>' % (
            '<br/>'.join(_tag['no_restric']), '<br/>'.join(_tag['no_multiple']), '<br/>'.join(_tag['stop']), '<br/>'.join(_tag['mut']), '<br/>'.join(good_ids), '<br/>'.join(_tag['amber'])))
    # Variable regions analysis
    w.write(
        '<div class="page-header"><a id="variable"></a><h2>Variable regions analysis</h2></div><p>The following group of sequences are identical clones on the variable '
        'regions:</p>')
    identical_clones_seq = _identical_clones.keys()
    if identical_clones_seq:
        for seq in identical_clones_seq:
            ids = list(set(_identical_clones[seq]))  # return only one occurrence of each item in the list
            w.write('<div class="row-fluid"><div class="span5"><pre>%d sequences (%.2f%% of valid sequences)<br/>%s</pre></div>' % (
                len(ids), float(len(ids)) / float(len(good_ids)) * 100, '<br/>'.join(ids)))
            w.write('<div class="span3"><table class="table table-striped table-bordered"><thead><tr><th>Variable region</th><th>Repeated sequence</th></tr></thead><tbody>')
            for z in range(len(_good_seq[ids[0]]['var'])):
                w.write('<td>%d</td><td>%s</td></tr>' % (z + 1, _good_seq[ids[0]]['var'][z]))
            w.write('</tbody></table></div></div>')
    else:
        w.write('<p>No clone was found.</p>')

    first = True
    for i in range(_nb_var_part):
        keys = []
        for k in _var_seq_common[str(i + 1)].keys():
            nb = _var_seq_common[str(i + 1)][k]
            if nb > 1:
                if first:
                    w.write(
                        '<p>Here\'s the distribution of the repeated sequences in variable regions:</p><table class="table table-striped table-bordered"><thead><tr><th>Variable '
                        'region</th><th>Repeated sequence</th><th>Number of occurrences (percentage of valid sequences)</th></tr></thead><tbody>')
                    first = False
                    keys.append(k)
                else:
                    keys.append(k)
        nb = len(keys)
        if nb != 0:
            w.write('<tr>')
            for z in range(nb):
                if z == 0:
                    w.write('<td rowspan="%d">%d</td>' % (nb, i + 1))
                w.write('<td>%s</td><td>%d (%.2f%%)</td></tr>' % (
                    keys[z], _var_seq_common[str(i + 1)][keys[z]], float(_var_seq_common[str(i + 1)][keys[z]]) / float(len(good_ids)) * 100))
    w.write('</tbody></table>')
    # Clustering
    w.write('<div class="page-header"><a id="cluster"></a><h2>Clustering</h2></div><p>The following clusters were generated by MCL:</p>')
    for line in open(mcl_output, 'r'):
        w.write('<div class="row-fluid"><div class="span6"><pre>%d sequences (%.2f%% of valid sequences)<br/>%s</pre></div></div>' % (
            len(line.split("\t")), float(len(line.split("\t"))) / float(len(good_ids)) * 100, '<br/>'.join(line.split("\t"))))
    # Statistics
    w.write('<div class="page-header"><a id="stat"></a><h2>Statistics</h2></div>')
    w.write('<p>Here\'s some statistics about the valid sequences:</p><p>Mean for the pairwise alignement scores: %.2f<br/>Standard deviation: %.2f</p>' % (
        float(numpy.mean(_align_scores)), float(numpy.std(_align_scores))))
    w.write('<div class="row-fluid"><div class="span6"><img src="%s" alt="Distribution of the pairwise alignment score"></div>' % get_os_path_name(graph_pic))
    w.write('<div class="span6"><table class="table table-striped table-bordered"><thead><tr><th>Pairwise Alignment Score</th><th>Number of occurrences</th></tr></thead><tbody>')
    uniq_scores = sorted(list(set(_align_scores)))
    scores_dic = {}
    for _score in uniq_scores:
        scores_dic[_score] = _align_scores.count(_score)
    scores_dic = dict(sorted(scores_dic.items()))
    scores = scores_dic.items()
    # scores.sort()
    for el in scores:
        w.write('<tr><td>%.2f</td><td>%d</td></tr>' % (el[0], el[1]))
    w.write('</tbody></table></div></div>')
    # Annex
    w.write('<div class="page-header"><a id="annex"></a><h2>Annex</h2></div>')
    w.write('<p><strong>Valid protein sequences</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
    for _id in good_ids:
        w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", _good_seq[_id]['prot'], re.DOTALL)))
    w.write('</textarea>')
    aln_out = generate_aln(_good_seq, good_ids)
    w.write(
        '<p>Multiple sequence alignment of the <strong>valid sequences</strong> generated by Clustal Omega:</p><textarea class="span8 fasta" type="text" rows="20" '
        'readonly="readonly">%s</textarea>' % str(
            aln_out))

    if _tag['no_multiple']:
        w.write(
            '<p><strong>Protein sequences with an incorrect number of nucleotides between the restriction sites</strong> in FASTA format:</p><textarea class="span8 fasta" '
            'type="text" rows="20" readonly="readonly">')
        for _id in _tag['no_multiple']:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", _all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')

    if _tag['mut']:
        w.write('<p><strong>Mutated protein sequences</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in _tag['mut']:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", _all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')
        aln_out = generate_aln(_all_seq_fasta, _tag['mut'])

        w.write(
            '<p>Multiple sequence alignment of the <strong>mutated sequences</strong> generated by Clustal Omega:</p><textarea class="span8 fasta" type="text" rows="20" '
            'readonly="readonly">%s</textarea>' % str(
                aln_out))

    if _tag['stop']:
        w.write('<p><strong>Protein sequences with a stop codon</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in _tag['stop']:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", _all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')

    if _tag['amber']:
        w.write('<p><strong>Protein sequences with an amber codon</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in _tag['amber']:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", _all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')

    w.write('</div></body></html>')
    w.close()


nb_seq = len(list(SeqIO.parse(args.input, "fasta")))

for seq_record in SeqIO.parse(args.input, "fasta"):
    seq_id = seq_record.id
    seq = str(seq_record.seq)
    seq = seq.upper()
    all_seq.append(seq_id)
    # Checking if both restriction sites are present in the sequence
    if site_res_5 in seq and site_res_3 in seq:
        valid = True
    else:
        valid = False
        tag['no_restric'].append(seq_id)
    # If sequence has both restriction sites, checking if it is necessary to take the reverse complement strand
    if valid:
        site_res_5_pos = seq.index(site_res_5)
        site_res_3_pos = seq.index(site_res_3)
        # If site_res_5_pos > site_res_3_pos, reverse complement strand has to be calculated
        if site_res_5_pos > site_res_3_pos:
            # Checking if the number of nucleic acids between the restriction sites is a multiple of 3
            length = math.fabs((site_res_5_pos + len(site_res_5)) - site_res_3_pos)
            valid = length % 3 == 0
            cut_seq = seq[:site_res_5_pos + len(site_res_5)]
            cut_seq = reverse_complement(cut_seq)

        # Else if site_res_5_pos < site_res_3_pos, use the sequence as it is
        else:
            # Checking if the number of nucleic acids between the restriction sites is a multiple of 3
            length = math.fabs((site_res_3_pos + len(site_res_3)) - site_res_5_pos)
            valid = length % 3 == 0
            cut_seq = seq[site_res_5_pos:]
        # If the number of nucleic acids between the restriction sites isn't a multiple of 3, put the sequence away
        if not valid:
            tag['no_multiple'].append(seq_id)
            prot_seq = translate(cut_seq)
            all_seq_fasta[seq_id] = {}
            all_seq_fasta[seq_id]['prot'] = prot_seq
        else:
            # Translate nucleic sequence into amino acid sequence
            prot_seq = translate(cut_seq)
            all_seq_fasta[seq_id] = {}
            all_seq_fasta[seq_id]['prot'] = prot_seq

            # Looking for stop codon in the sequence and getting their position in the sequence
            if '*' in prot_seq:
                pos_stop = [m.start() for m in re.finditer(r"\*", prot_seq)]
                stop = False
                # Checking if stop codon is between the restriction sites, also checking if it is an amber codon. if stop codon other than amber codon -> tag stop
                for i in range(len(pos_stop)):
                    if pos_stop[i] < length / 3:
                        stop_codon_nuc = cut_seq[pos_stop[i] * 3:pos_stop[i] * 3 + 3]
                        if stop_codon_nuc != "TAG":
                            tag['stop'].append(seq_id)
                            stop = True
                            break
                        else:
                            if seq_id not in tag['amber']:
                                tag['amber'].append(seq_id)
                # If stop codon wasn't found between the restriction sites
                if not stop:
                    """
                    # Checking if there is a stop codon outside the restriction sites. If yes -> tag ok_stop_ext
                    for i in range(len(pos_stop)):
                        if (pos_stop[i] > length/3):
                            stop_codon_nuc = cut_seq[pos_stop[i]*3:pos_stop[i]*3+3]
                            if stop_codon_nuc != "TAG":
                                tag['ok_stop_ext'].append(seq_id)
                                stop = True
                                break
                            else:
                                if (seq_id not in tag['amber']):
                                    tag['amber'].append(seq_id)
                    """
                    # Checking if there was a mutation in the fix part, if yes -> tag mut else retrieve variable parts
                    mut = False
                    pattern_part = args.pattern.split(":")
                    tmp_prot_seq = prot_seq
                    var_parts = []
                    for i in range(len(pattern_part) - 1):  # not checking the latest fix part
                        part = pattern_part[i]
                        # If part is fix
                        if not part[0].isdigit():
                            # If part not in prot_seq -> mutation, flag then break
                            if part not in tmp_prot_seq:
                                mut = True
                                tag['mut'].append(seq_id)
                                break
                            # Else, store the variable part if exist then remove the fix part + variable part (tmp_prot_seq starts at the end of part)
                            else:
                                pos_fix = tmp_prot_seq.index(part)
                                if pos_fix != 0:
                                    var_parts.append(tmp_prot_seq[0:pos_fix])
                                tmp_prot_seq = tmp_prot_seq[pos_fix + len(part):]
                        # Else part is variable
                        else:
                            nb_var_part += 1
                    # Treating latest fix part if no mutation before
                    if not mut:
                        last_part = pattern_part[-1]
                        last_var = pattern_part[-2]
                        if '-' in last_var:
                            var_max = int(last_var.split('-')[1])
                        else:
                            var_max = int(last_var)
                        last_part = last_part[0:var_max + 1]
                        if last_part not in tmp_prot_seq:
                            mut = True
                            tag['mut'].append(seq_id)
                        else:
                            pos_fix = tmp_prot_seq.index(last_part)
                            if pos_fix != 0:
                                var_parts.append(tmp_prot_seq[0:pos_fix])
                    # If no mutation the sequence is validated and all the info are stored
                    if not mut:
                        good_seq[seq_id] = {}
                        good_seq[seq_id]['dna'] = cut_seq
                        good_seq[seq_id]['prot'] = prot_seq
                        good_seq[seq_id]['var'] = var_parts

# If all sequences are invalid, the program will exit as there is no data to continue
if not good_seq:
    sys.exit("There is only one valid sequence among the input data. At least 2 valid sequences are necessary to proceed to the next step. The program will now exit")
elif len(good_seq.keys()) == 1:

    sys.exit("There is only one valid sequence among the input data. At least 2 valid sequences are necessary to proceed to the next step. The program will now exit")

# Initialization of dict var_seq_common
for n in range(nb_var_part):
    var_seq_common[str(n + 1)] = {}

# Opening the file where the mcl input will be written
with open(mcl_file, 'w+') as mcl:
    seq_keys = good_seq.keys()
    for i in range(len(seq_keys)):
        var_1 = good_seq[list(seq_keys)[i]]['var']

        # Classifying variable sequences
        for k in range(len(var_1)):
            try:
                var_seq_common[str(k + 1)][var_1[k]] += 1
            except KeyError:
                var_seq_common[str(k + 1)][var_1[k]] = 1

        for j in range(i + 1, len(seq_keys)):
            var_2 = good_seq[list(seq_keys)[j]]['var']
            score = 0.0
            # Comparing the sequences' variable parts to find identical clones
            if var_1 == var_2:
                try:
                    clone_seq = "".join(var_1)
                    identical_clones[clone_seq].extend([seq_keys[i], seq_keys[j]])
                except KeyError:
                    identical_clones[clone_seq] = [seq_keys[i], seq_keys[j]]
            # Align the 2 sequences using NWalign_PAM30 => replace by pairwise2
            seq_1 = ''.join(var_1)
            seq_2 = ''.join(var_2)
            matrix = MatrixInfo.pam30
            if len(seq_2) > len(seq_1):
                score = get_identity(pairwise2.align.globalds(seq_1, seq_2, matrix, -11, -1)[0][0], pairwise2.align.globalds(seq_1, seq_2, matrix, -11, -1)[0][1]) * 100
            else:
                score = get_identity(pairwise2.align.globalds(seq_2, seq_1, matrix, -11, -1)[0][0], pairwise2.align.globalds(seq_2, seq_1, matrix, -11, -1)[0][1]) * 100
            align_scores.append(score)
            mcl.write('%s\t%s\t%0.2f\n' % (list(seq_keys)[i], list(seq_keys)[j], score))

# Clusters formation
subprocess.call(["mcl", mcl_file, "--abc", "-I", "6.0", "-o", mcl_output], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Producing distribution graph
plot.hist(align_scores, bins=numpy.arange(0, 101, 2))
plot.xlabel('Pairwise Alignment Score')
plot.ylabel('Number of occurrences')
plot.title('Distribution of the pairwise alignment score')
plot.grid(True)
plot.savefig(graph_pic)

# Generating html report
report_html(html_file, tag, all_seq, good_seq, all_seq_fasta, identical_clones, nb_var_part, var_seq_common, align_scores, args)

# Removing intermediate files
subprocess.call(["rm", mcl_file, mcl_output], shell=False)
print("HTML report has been generated in the output directory. The program will now exit.")
