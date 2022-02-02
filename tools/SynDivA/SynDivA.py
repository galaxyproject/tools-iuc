# !/usr/bin/env python title           :SynDivA.py description     :This script will analyze fasta files, look for restriction sites, cut the sequences around the restriction
# sites, translate the nucleic sequences into amino acids sequences. author          :Fabienne Wong Jun Tai date            :20121107 version         :1.0 usage
# :python SynDivA.py -i file.fasta -o /output/dir/ -p pattern -5 seq_restric_5'-3 seq_restric_3' notes           : python_version  :3.7.11 biopython_max_version  :1.72
# ==============================================================================
import re
import subprocess

import math
import matplotlib
import matplotlib.pyplot as plot
import numpy
from Bio import SeqIO, Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment

from args import *

matplotlib.use('Agg')

args = Args()
print(sys.path[0])
# Variables initialization
SynDivA_script_dir = sys.path[0]
print(SynDivA_script_dir)
directory = args.output_dir
mcl_file = directory + "mcl.in"
mcl_output = directory + "mcl.out"
html_file = directory + "SynDivA_report.html"
graph_pic = directory + "distri.png"
input_file = os.path.basename(args.input)
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


def reverse_complement(seq):
    # Generate the reverse complement
    complement_nuc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_com = ""
    for n in (seq[::-1]):
        rev_com += complement_nuc[n]
    return rev_com


def generate_aln(seq_dic, ids):
    # Multiple Sequence Alignment via ClustalO
    input = ''
    for k in ids:
        input += '>%s\n%s\n' % (k, re.sub("(.{80})", "\\1\n", seq_dic[k]['prot'], re.DOTALL))
    p = subprocess.Popen("clustalo -i - --outfmt clu", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
    aln_out, aln_err = p.communicate(input=input)
    print(type(aln_out))
    return aln_out


def report_html(html_file, tag, all_seq, good_seq, all_seq_fasta, identical_clones, nb_var_part, var_seq_common, align_scores, args):
    # Generate the html file for the report
    all_seq.sort()
    no_restric = tag['no_restric']
    no_restric.sort()
    no_multiple = tag['no_multiple']
    no_multiple.sort()
    stop = tag['stop']
    stop.sort()
    amber = tag['amber']
    amber.sort()
    mut = tag['mut']
    mut.sort()
    # good_ids = good_seq.keys()

    good_seq = dict(sorted(good_seq.items()))
    good_ids = good_seq.keys()

    # good_ids.sort()

    w = open(html_file, 'w')
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
            input_file, len(all_seq), args.pattern, args.site_res_5, args.site_res_3))
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
            len(no_restric), float(len(no_restric)) / float(len(all_seq)) * 100, len(no_multiple), float(len(no_multiple)) / float(len(all_seq)) * 100, len(stop),
            float(len(stop)) / float(len(all_seq)) * 100, len(mut), float(len(mut)) / float(len(all_seq)) * 100, len(good_ids), float(len(good_ids)) / float(len(all_seq)) * 100,
            len(amber)))
    w.write(
        '<tr><td class="text-error">%s</td><td class="text-error">%s</td><td class="text-error">%s</td><td class="text-warning">%s</td><td '
        'class="text-success">%s</td><td>%s</td></tr></table>' % (
            '<br/>'.join(no_restric), '<br/>'.join(no_multiple), '<br/>'.join(stop), '<br/>'.join(mut), '<br/>'.join(good_ids), '<br/>'.join(amber)))
    # Variable regions analysis
    w.write(
        '<div class="page-header"><a id="variable"></a><h2>Variable regions analysis</h2></div><p>The following group of sequences are identical clones on the variable '
        'regions:</p>')
    identical_clones_seq = identical_clones.keys()
    if identical_clones_seq:
        for seq in identical_clones_seq:
            ids = list(set(identical_clones[seq]))  # return only one occurrence of each item in the list
            w.write('<div class="row-fluid"><div class="span5"><pre>%d sequences (%.2f%% of valid sequences)<br/>%s</pre></div>' % (
                len(ids), float(len(ids)) / float(len(good_ids)) * 100, '<br/>'.join(ids)))
            w.write('<div class="span3"><table class="table table-striped table-bordered"><thead><tr><th>Variable region</th><th>Repeated sequence</th></tr></thead><tbody>')
            for z in range(len(good_seq[ids[0]]['var'])):
                w.write('<td>%d</td><td>%s</td></tr>' % (z + 1, good_seq[ids[0]]['var'][z]))
            w.write('</tbody></table></div></div>')
    else:
        w.write('<p>No clone was found.</p>')

    first = True
    for i in range(nb_var_part):
        keys = []
        for k in (var_seq_common[str(i + 1)].keys()):
            nb = var_seq_common[str(i + 1)][k]
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
                    keys[z], var_seq_common[str(i + 1)][keys[z]], float(var_seq_common[str(i + 1)][keys[z]]) / float(len(good_ids)) * 100))
    w.write('</tbody></table>')
    # Clustering
    w.write('<div class="page-header"><a id="cluster"></a><h2>Clustering</h2></div><p>The following clusters were generated by MCL:</p>')
    for line in open(mcl_output, 'r'):
        w.write('<div class="row-fluid"><div class="span6"><pre>%d sequences (%.2f%% of valid sequences)<br/>%s</pre></div></div>' % (
            len(line.split("\t")), float(len(line.split("\t"))) / float(len(good_ids)) * 100, '<br/>'.join(line.split("\t"))))
    # Statistics
    w.write('<div class="page-header"><a id="stat"></a><h2>Statistics</h2></div>')
    w.write('<p>Here\'s some statistics about the valid sequences:</p><p>Mean for the pairwise alignement scores: %.2f<br/>Standard deviation: %.2f</p>' % (
        numpy.mean(align_scores), numpy.std(align_scores)))
    w.write('<div class="row-fluid"><div class="span6"><img src="%s" alt="Distribution of the pairwise alignment score"></div>' % os.path.basename(graph_pic))
    w.write('<div class="span6"><table class="table table-striped table-bordered"><thead><tr><th>Pairwise Alignment Score</th><th>Number of occurrences</th></tr></thead><tbody>')
    uniq_scores = sorted(list(set(align_scores)))
    scores_dic = {}
    for score in uniq_scores:
        scores_dic[score] = align_scores.count(score)

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
        w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", good_seq[_id]['prot'], re.DOTALL)))
    w.write('</textarea>')
    aln_out = generate_aln(good_seq, good_ids)
    print(str(aln_out))
    w.write(
        '<p>Multiple sequence alignment of the <strong>valid sequences</strong> generated by Clustal Omega:</p><textarea class="span8 fasta" type="text" rows="20" '
        'readonly="readonly">%s</textarea>' % str(
            aln_out))

    if no_multiple:
        w.write(
            '<p><strong>Protein sequences with an incorrect number of nucleotides between the restriction sites</strong> in FASTA format:</p><textarea class="span8 fasta" '
            'type="text" rows="20" readonly="readonly">')
        for _id in no_multiple:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')

    if mut:
        w.write('<p><strong>Mutated protein sequences</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in mut:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')
        aln_out = generate_aln(all_seq_fasta, mut)

        w.write(
            '<p>Multiple sequence alignment of the <strong>mutated sequences</strong> generated by Clustal Omega:</p><textarea class="span8 fasta" type="text" rows="20" '
            'readonly="readonly">%s</textarea>' % str(
                aln_out))

    if stop:
        w.write('<p><strong>Protein sequences with a stop codon</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in stop:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", all_seq_fasta[_id]['prot'], re.DOTALL)))
        w.write('</textarea>')

    if amber:
        w.write('<p><strong>Protein sequences with an amber codon</strong> in FASTA format:</p><textarea class="span8 fasta" type="text" rows="20" readonly="readonly">')
        for _id in amber:
            w.write('>%s\n%s\n' % (_id, re.sub("(.{80})", "\\1\n", all_seq_fasta[_id]['prot'], re.DOTALL)))
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
            prot_seq = Seq.translate(cut_seq)
            all_seq_fasta[seq_id] = {}
            all_seq_fasta[seq_id]['prot'] = prot_seq
        else:
            # Translate nucleic sequence into amino acid sequence
            prot_seq = Seq.translate(cut_seq)
            all_seq_fasta[seq_id] = {}
            all_seq_fasta[seq_id]['prot'] = prot_seq

            # Looking for stop codon in the sequence and getting their position in the sequence
            if '*' in prot_seq:
                pos_stop = [m.start() for m in re.finditer("\*", prot_seq)]
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
    print("All sequences are invalid. At least 2 valid sequences are necessary to proceed to the next step. The program will now exit.")
    sys.exit()
elif len(good_seq.keys()) == 1:
    print("There is only one valid sequence among the input data. At least 2 valid sequences are necessary to proceed to the next step. The program will now exit")
    sys.exit()

# Initialization of dict var_seq_common
for n in range(nb_var_part):
    var_seq_common[str(n + 1)] = {}

# Opening the file where the mcl input will be written
mcl = open(mcl_file, 'w')

id = good_seq.keys()
for i in range(len(id)):
    var_1 = good_seq[list(id)[i]]['var']

    # Classifying variable sequences
    for k in range(len(var_1)):
        try:
            var_seq_common[str(k + 1)][var_1[k]] += 1
        except KeyError:
            var_seq_common[str(k + 1)][var_1[k]] = 1

    for j in range(i + 1, len(id)):
        var_2 = good_seq[list(id)[j]]['var']
        # Comparing the sequences' variable parts to find identical clones
        if var_1 == var_2:
            try:
                s = "".join(var_1)
                identical_clones[s].extend([id[i], id[j]])
            except KeyError:
                identical_clones[s] = [id[i], id[j]]

        # Align the 2 sequences using NWalign_PAM30
        seq_1 = ''.join(var_1)
        seq_2 = ''.join(var_2)
        print(seq_1)
        print(seq_2)
        matrix = matlist.pam30
        cpt = 0
        if len(seq_2) > len(seq_1):
            print(pairwise2.align.globalds(seq_1, seq_2, matrix, -11, -1))
            for a in pairwise2.align.globalds(seq_1, seq_2, matrix, -11, -1):
                for k in range(a[4]):
                    if a[0][k] == a[1][k]:
                        cpt += 1
                print(format_alignment(*a, full_sequences=True))
        else:
            print(pairwise2.align.globalds(seq_2, seq_1, matrix, -11, -1))
            for a in pairwise2.align.globalds(seq_2, seq_1, matrix, -11, -1):
                for k in range(a[4]):
                    if a[0][k] == a[1][k]:
                        cpt += 1
                print(format_alignment(*a, full_sequences=True))
        print("######################################@")
        print(cpt)

        if len(seq_2) > len(seq_1):
            p = subprocess.Popen(SynDivA_script_dir + "/NWalign_PAM30 %s %s 3" % (seq_1, seq_2), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen(SynDivA_script_dir + "/NWalign_PAM30 %s %s 3" % (seq_2, seq_1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        out, err = p.communicate()

        print(out)
        print("######################################@")
        lines = out.split(bytes("\n", encoding='utf8'))
        print(lines[5].split(bytes(' ', encoding='utf8'))[5])
        score = float(lines[5].split(bytes(' ', encoding='utf8'))[5]) * 100
        align_scores.append(score)
        mcl.write('%s\t%s\t%0.2f\n' % (list(id)[i], list(id)[j], score))
mcl.close()

# Clusters formation
subprocess.call("mcl %s --abc -I 6.0 -o %s" % (mcl_file, mcl_output), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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
subprocess.call("rm %s %s " % (mcl_file, mcl_output), shell=True)

print("HTML report has been generated in the output directory. The program will now exit.")
