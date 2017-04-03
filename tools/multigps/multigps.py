import argparse
import os
import shutil
import subprocess
import sys
import tempfile

BUFF_SIZE = 1048576
DESIGN_FILE = 'design.tabular'

parser = argparse.ArgumentParser()
parser.add_argument('--all_events_table', dest='all_events_table', help='Output all events table file')
parser.add_argument('--alphascale', dest='alphascale', type=float, default=None, help='Alpha scaling factor')
parser.add_argument('--chrom_len_file', dest='chrom_len_file', help='File listing the lengths of all chromosomes')
parser.add_argument('--control', dest='control', default=None, help='Input control files and data formats')
parser.add_argument('--diffp', dest='diffp', type=float, default=None, help='Minimum p-value for reporting differential enrichment')
parser.add_argument('--edgerod', dest='edgerod', type=float, default=None, help='EdgeR over-dispersion parameter value')
parser.add_argument('--exclude', dest='exclude', default=None, help='File containing a set of regions to ignore during MultiGPS training')
parser.add_argument('--expt', dest='expt', default=None, help="Input expt files and data formats")
parser.add_argument('--eventsaretxt', dest='eventsaretxt', default=None, help='Append a .txt extension to the events file for browser rendering')
parser.add_argument('--fixedalpha', dest='fixedalpha', type=int, default=None, help='Impose this alpha')
parser.add_argument('--fixedmodelrange', dest='fixedmodelrange', default=None, help='Keep binding model range fixed to inital size')
parser.add_argument('--fixedpb', dest='fixedpb', type=int, default=None, help='Fixed per-base limit')
parser.add_argument('--fixedscaling', dest='fixedscaling', type=float, default=None, help='Multiply control counts by total tag count ratio and then by this factor')
parser.add_argument('--format', dest='format', default=None, help='Input expt file data format')
parser.add_argument('--gaussmodelsmoothing', dest='gaussmodelsmoothing', default=None, help='Use Gaussian model smoothing')
parser.add_argument('--gausssmoothparam', dest='gausssmoothparam', type=int, default=None, help='Smoothing factor')
parser.add_argument('--input_item', dest='input_items', action='append', nargs=7, default=None, help="Input files, attributes and options")
parser.add_argument('--jointinmodel', dest='jointinmodel', default=None, help='Allow joint events in model updates')
parser.add_argument('--mappability', dest='mappability', type=float, default=None, help='Fraction of the genome that is mappable for these experiments')
parser.add_argument('--maxtrainingrounds', dest='maxtrainingrounds', type=int, default=None, help='Maximum number of training rounds for updating binding event read distributions')
parser.add_argument('--medianscale', dest='medianscale', default=None, help='Use the median signal/control ratio as the scaling factor')
parser.add_argument('--meme1proc', dest='meme1proc', default=None, help='Do not run the parallel version of meme')
parser.add_argument('--mememaxw', dest='mememaxw', type=int, default=None, help='Maximum motif width for MEME')
parser.add_argument('--mememinw', dest='mememinw', type=int, default=None, help='Minimum motif width for MEME')
parser.add_argument('--memenmotifs', dest='memenmotifs', type=int, default=None, help='Number of motifs MEME should find for each condition')
parser.add_argument('--minfold', dest='minfold', type=float, default=None, help='Minimum event fold-change vs scaled control')
parser.add_argument('--minqvalue', dest='minqvalue', type=float, default=None, help='Minimum Q-value (corrected p-value) of reported binding events')
parser.add_argument('--minmodelupdateevents', dest='minmodelupdateevents', type=int, default=None, help='Minimum number of events to support an update of the read distribution')
parser.add_argument('--mlconfignotshared', dest='mlconfignotshared', default=None, help='Share component configs in the ML step')
parser.add_argument('--nocache', dest='nocache', default=None, help='Turn off caching of the entire set of experiments')
parser.add_argument('--nodifftests', dest='nodifftests', default=None, help='Run differential enrichment tests')
parser.add_argument('--nomodelsmoothing', dest='nomodelsmoothing', default=None, help='Perform binding model smoothing')
parser.add_argument('--nomodelupdate', dest='nomodelupdate', default=None, help='Perform binding model updates')
parser.add_argument('--nomotifprior', dest='nomotifprior', default=None, help='Perform motif-finding only')
parser.add_argument('--nomotifs', dest='nomotifs', default=None, help='Perform motif-finding and motif priors')
parser.add_argument('--nonunique', dest='nonunique', default=None, help='Use non-unique reads')
parser.add_argument('--noposprior', dest='noposprior', default=None, help='Perform inter-experiment positional prior')
parser.add_argument('--noscaling', dest='noscaling', default=None, help='Do not use signal vs control scaling')
parser.add_argument('--output_bed', dest='output_bed', help='Output bed results file')
parser.add_argument('--output_html', dest='output_html', help='Output html results file')
parser.add_argument('--output_html_files_path', dest='output_html_files_path', help='Output html extra files')
parser.add_argument('--plotscaling', dest='plotscaling', default=None, help='Plot diagnostic information for the chosen scaling method')
parser.add_argument('--poissongausspb', dest='poissongausspb', type=int, default=None, help='Poisson threshold for filtering per base')
parser.add_argument('--prlogconf', dest='prlogconf', type=int, default=None, help='Poisson log threshold for potential region scanning')
parser.add_argument('--probshared', dest='probshared', type=float, default=None, help='Probability that events are shared across conditions')
parser.add_argument('--readdistributionfile', dest='readdistributionfile', default=None, help='Optional binding event read distribution file for initializing models')
parser.add_argument('--regressionscale', dest='regressionscale', default=None, help='Use scaling by regression on binned tag counts')
parser.add_argument('--replicates_counts', dest='replicates_counts', help='Output replicates counts file')
parser.add_argument('--scalewin', dest='scalewin', type=int, default=None, help='Window size for estimating scaling ratios')
parser.add_argument('--seq', dest='seq', default=None, help='Reference genome path')
parser.add_argument('--sesscale', dest='sesscale', default=None, help='Estimate scaling factor by SES')
parser.add_argument('--splinesmoothparam', dest='splinesmoothparam', type=int, default=None, help='Spline smoothing parameter')
parser.add_argument('--threads', dest='threads', type=int, default=4, help='The number of threads to run')
args = parser.parse_args()


def generate_design_file(input_items):
    design_file = open(DESIGN_FILE, 'w')
    for item in input_items:
        file_name, label, file_format, condition_name, replicate_name, experiment_type, fixed_read_count = item
        file_format = file_format.upper()
        items = [file_name, label, file_format, condition_name]
        if replicate_name not in ['None', None, '']:
            items.append(replicate_name)
        if experiment_type not in ['None', None, '']:
            items.append(experiment_type)
        if fixed_read_count not in ['None', None, '']:
            items.append(fixed_read_count)
        design_file.write('%s\n' % '\t'.join(items))
    design_file.close()


def get_file_with_extension(dir, ext):
    file_list = [f for f in os.listdir(dir) if f.endswith(ext)]
    if len(file_list) == 1:
        return file_list[0]
    stop_err('Error running MultiGPS: output file with extension "%s" not generated.' % ext)


def get_stderr_exception(tmp_err, tmp_stderr, tmp_out, tmp_stdout, include_stdout=False):
    tmp_stderr.close()
    # Get stderr, allowing for case where it's very large.
    tmp_stderr = open(tmp_err, 'rb')
    stderr_str = ''
    buffsize = BUFF_SIZE
    try:
        while True:
            stderr_str += tmp_stderr.read(buffsize)
            if not stderr_str or len(stderr_str) % buffsize != 0:
                break
    except OverflowError:
        pass
    tmp_stderr.close()
    if include_stdout:
        tmp_stdout = open(tmp_out, 'rb')
        stdout_str = ''
        buffsize = BUFF_SIZE
        try:
            while True:
                stdout_str += tmp_stdout.read(buffsize)
                if not stdout_str or len(stdout_str) % buffsize != 0:
                    break
        except OverflowError:
            pass
    tmp_stdout.close()
    if include_stdout:
        return 'STDOUT\n%s\n\nSTDERR\n%s\n' % (stdout_str, stderr_str)
    return stderr_str


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


# Preparation.
tmp_dir = tempfile.mkdtemp(prefix="tmp-multigps-")
os.makedirs(args.output_html_files_path)
# Build the command line.
cmd = 'multigps'
# General options
cmd += ' --threads %s' % args.threads
if args.eventsaretxt is not None:
    # Append .txt extensions to events hrefs
    # in output dataset so files will render
    # in the browser.
    cmd += ' --eventsaretxt'
if args.meme1proc is not None:
    # Do not run the parallel version of meme.
    cmd += ' --meme1proc'
# Experiment.
if args.input_items is not None:
    generate_design_file(args.input_items)
    cmd += ' --design %s' % DESIGN_FILE
else:
    if args.expt is not None:
        cmd += ' --expt %s' % args.expt
    if args.format is not None:
        cmd += ' --format %s' % args.format
    if args.control is not None:
        cmd += ' --ctrl %s' % args.control
cmd += ' --geninfo %s' % args.chrom_len_file
# Advanced options.
if args.seq is not None:
    cmd += ' --seq %s' % args.seq
# Limits on how many reads
if args.fixedpb is not None:
    cmd += ' --fixedpb %d' % args.fixedpb
if args.poissongausspb is not None:
    cmd += ' --poissongausspb %d' % args.poissongausspb
if args.nonunique is not None:
    cmd += ' --nonunique'
if args.mappability is not None:
    cmd += ' --mappability %4f' % args.mappability
if args.nocache is not None:
    cmd += ' --nocache'
# Scaling data.'
if args.noscaling is not None:
    cmd += ' --noscaling %s' % args.noscaling
if args.medianscale is not None:
    cmd += ' --medianscale %s' % args.medianscale
if args.regressionscale is not None:
    cmd += ' --regressionscale %s' % args.regressionscale
if args.sesscale is not None:
    cmd += ' --sesscale %s' % args.sesscale
if args.fixedscaling is not None:
    cmd += ' --fixedscaling %4f' % args.fixedscaling
if args.scalewin is not None:
    cmd += ' --scalewin %d' % args.scalewin
if args.plotscaling is not None:
    cmd += ' --plotscaling %s' % args.plotscaling
# Running MultiGPS.
if args.readdistributionfile is not None:
    cmd += ' --d %s' % args.readdistributionfile
if args.maxtrainingrounds is not None:
    cmd += ' --r %s' % args.maxtrainingrounds
if args.nomodelupdate is not None:
    cmd += ' --nomodelupdate'
if args.minmodelupdateevents is not None:
    cmd += ' --minmodelupdateevents %d' % args.minmodelupdateevents
if args.nomodelsmoothing is not None:
    cmd += ' --nomodelsmoothing'
if args.splinesmoothparam is not None:
    cmd += ' --splinesmoothparam %d' % args.splinesmoothparam
if args.gaussmodelsmoothing is not None:
    cmd += ' --gaussmodelsmoothing'
if args.gausssmoothparam is not None:
    cmd += ' --gausssmoothparam %s' % args.gausssmoothparam
if args.jointinmodel is not None:
    cmd += ' --jointinmodel'
if args.fixedmodelrange is not None:
    cmd += ' --fixedmodelrange'
if args.prlogconf is not None:
    cmd += ' --prlogconf %d' % args.prlogconf
if args.fixedalpha is not None:
    cmd += ' --fixedalpha %d' % args.fixedalpha
if args.alphascale is not None:
    cmd += ' --alphascale %4f' % args.alphascale
if args.mlconfignotshared is not None:
    cmd += ' --mlconfignotshared'
if args.exclude not in [None, 'None']:
    cmd += ' --exclude %s' % args.exclude_file
# MultiGPS priors
if args.noposprior is not None:
    cmd += ' --noposprior'
if args.probshared is not None:
    cmd += ' --probshared %4f' % args.probshared
if args.memenmotifs is not None:
    cmd += ' --memenmotifs %d' % args.memenmotifs
if args.mememinw is not None:
    cmd += ' --mememinw %d' % args.mememinw
if args.mememaxw is not None:
    cmd += ' --mememaxw %d' % args.mememaxw
if args.nomotifs is not None:
    cmd += ' --nomotifs'
if args.nomotifprior is not None:
    cmd += ' --nomotifprior'
# Reporting binding events
if args.minqvalue is not None:
    cmd += ' --q %4f' % args.minqvalue
if args.minfold is not None:
    cmd += ' --minfold %4f' % args.minfold
if args.nodifftests is not None:
    cmd += ' --nodifftests'
if args.edgerod is not None:
    cmd += ' --edgerod %4f' % args.edgerod
if args.diffp is not None:
    cmd += ' --diffp %4f' % args.diffp
# Output directory.
cmd += ' --out %s' % args.output_html_files_path
# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile(dir=tmp_dir).name
tmp_stdout = open(tmp_out, 'wb')
tmp_err = tempfile.NamedTemporaryFile(dir=tmp_dir).name
tmp_stderr = open(tmp_err, 'wb')
tmp_filename = tempfile.NamedTemporaryFile(dir=tmp_dir).name
# Execute the command.
proc = subprocess.Popen(args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True)
rc = proc.wait()
if rc != 0:
    error_message = get_stderr_exception(tmp_err, tmp_stderr, tmp_out, tmp_stdout)
    stop_err(error_message)
# Move each output file to the approapriate output dataset path.
output_bed = get_file_with_extension(args.output_html_files_path, 'bed')
shutil.move(os.path.join(args.output_html_files_path, output_bed), args.output_bed)
output_html = get_file_with_extension(args.output_html_files_path, 'html')
shutil.move(os.path.join(args.output_html_files_path, output_html), args.output_html)
replicates_counts = get_file_with_extension(args.output_html_files_path, 'counts')
shutil.move(os.path.join(args.output_html_files_path, replicates_counts), args.replicates_counts)
all_events_table = get_file_with_extension(args.output_html_files_path, 'table.txt')
shutil.move(os.path.join(args.output_html_files_path, all_events_table), args.all_events_table)
# Clean up.
if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
