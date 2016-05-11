#!/usr/bin/env python
import argparse
import os
import shutil
import string
import subprocess
import sys
import tempfile

BUFFSIZE = 1048576
# Translation table for reverse Complement, with ambiguity codes.
DNA_COMPLEMENT = string.maketrans("ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb")


def reverse(sequence):
    # Reverse sequence string.
    return sequence[::-1]


def dna_complement(sequence):
    # Complement DNA sequence string.
    return sequence.translate(DNA_COMPLEMENT)


def dna_reverse_complement(sequence):
    # Returns the reverse complement of the sequence.
    sequence = reverse(sequence)
    return dna_complement(sequence)


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument('--input_motifs', dest='input_motifs', help='MEME output formatted files for input to fimo')
parser.add_argument('--input_fasta', dest='input_fasta', help='Fassta sequence file')
parser.add_argument('--options_type', dest='options_type', help='Basic or Advance options')
parser.add_argument('--input_psp', dest='input_psp', default=None, help='File containing position specific priors')
parser.add_argument('--input_prior_dist', dest='input_prior_dist', default=None, help='File containing binned distribution of priors')
parser.add_argument('--alpha', dest='alpha', type=float, default=1.0, help='The alpha parameter for calculating position specific priors')
parser.add_argument('--bgfile', dest='bgfile', default=None, help='Background file type, used only if not "default"')
parser.add_argument('--max_strand', action='store_true', help='If matches on both strands at a given position satisfy the output threshold, only report the match for the strand with the higher score')
parser.add_argument('--max_stored_scores', dest='max_stored_scores', type=int, help='Maximum score count to store')
parser.add_argument('--motif', dest='motifs', action='append', default=[], help='Specify motif by id')
parser.add_argument('--motif_pseudo', dest='motif_pseudo', type=float, default=0.1, help='Pseudocount to add to counts in motif matrix')
parser.add_argument('--no_qvalue', action='store_true', help='Do not compute a q-value for each p-value')
parser.add_argument('--norc', action='store_true', help='Do not score the reverse complement DNA strand')
parser.add_argument('--output_path', dest='output_path', help='Output files directory')
parser.add_argument('--parse_genomic_coord', action='store_true', help='Check each sequence header for UCSC style genomic coordinates')
parser.add_argument('--qv_thresh', action='store_true', help='Use q-values for the output threshold')
parser.add_argument('--thresh', dest='thresh', type=float, help='p-value threshold')
parser.add_argument('--gff_output', dest='gff_output', help='Gff output file')
parser.add_argument('--html_output', dest='html_output', help='HTML output file')
parser.add_argument('--interval_output', dest='interval_output', help='Interval output file')
parser.add_argument('--txt_output', dest='txt_output', help='Text output file')
parser.add_argument('--xml_output', dest='xml_output', help='XML output file')
args = parser.parse_args()

fimo_cmd_list = ['fimo']
if args.options_type == 'advanced':
    fimo_cmd_list.append('--alpha %4f' % args.alpha)
    if args.bgfile is not None:
        fimo_cmd_list.append('--bgfile "%s"' % args.bgfile)
    if args.max_strand:
        fimo_cmd_list.append('--max-strand')
    fimo_cmd_list.append('--max-stored-scores %d' % args.max_stored_scores)
    if len(args.motifs) > 0:
        for motif in args.motifs:
            fimo_cmd_list.append('--motif "%s"' % motif)
    fimo_cmd_list.append('--motif-pseudo %4f' % args.motif_pseudo)
    if args.no_qvalue:
        fimo_cmd_list.append('--no-qvalue')
    if args.norc:
        fimo_cmd_list.append('--norc')
    if args.parse_genomic_coord:
        fimo_cmd_list.append('--parse-genomic-coord')
    if args.qv_thresh:
        fimo_cmd_list.append('--qv-thresh')
    fimo_cmd_list.append('--thresh %4f' % args.thresh)
    if args.input_psp is not None:
        fimo_cmd_list.append('--psp "%s"' % args.input_psp)
    if args.input_prior_dist is not None:
        fimo_cmd_list.append('--prior-dist "%s"' % args.input_prior_dist)
fimo_cmd_list.append('--o "%s"' % (args.output_path))
fimo_cmd_list.append('--verbosity 1')
fimo_cmd_list.append(args.input_motifs)
fimo_cmd_list.append(args.input_fasta)

fimo_cmd = ' '.join(fimo_cmd_list)

try:
    tmp_stderr = tempfile.NamedTemporaryFile()
    proc = subprocess.Popen(args=fimo_cmd, shell=True, stderr=tmp_stderr)
    returncode = proc.wait()
    tmp_stderr.seek(0)
    stderr = ''
    try:
        while True:
            stderr += tmp_stderr.read(BUFFSIZE)
            if not stderr or len(stderr) % BUFFSIZE != 0:
                break
    except OverflowError:
        pass
    if returncode != 0:
        stop_err(stderr)
except Exception as e:
    stop_err('Error running FIMO:\n%s' % str(e))

shutil.move(os.path.join(args.output_path, 'fimo.txt'), args.txt_output)
shutil.move(os.path.join(args.output_path, 'fimo.gff'), args.gff_output)
shutil.move(os.path.join(args.output_path, 'fimo.xml'), args.xml_output)
shutil.move(os.path.join(args.output_path, 'fimo.html'), args.html_output)

out_file = open(args.interval_output, 'wb')
out_file.write("#%s\n" % "\t".join(("chr", "start", "end", "pattern name", "score", "strand", "matched sequence", "p-value", "q-value")))
for line in open(args.txt_output):
    if line.startswith('#'):
        continue
    fields = line.rstrip("\n\r").split("\t")
    start, end = int(fields[2]), int(fields[3])
    sequence = fields[7]
    if start > end:
        # Flip start and end and set strand.
        start, end = end, start
        strand = "-"
        # We want sequences relative to strand; FIMO always provides + stranded sequence.
        sequence = dna_reverse_complement(sequence)
    else:
        strand = "+"
    # Make 0-based start position.
    start -= 1
    out_file.write("%s\n" % "\t".join([fields[1], str(start), str(end), fields[0], fields[4], strand, sequence, fields[5], fields[6]]))
out_file.close()
