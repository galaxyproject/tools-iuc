#!/usr/bin/env python

# Gets interfaced by Galaxy or can be used for bash scripting
import argparse
import logging
import os
import sys

import output_report
import spyboat
from numpy import float32
from skimage import io

logging.basicConfig(level=logging.INFO, stream=sys.stdout, force=True)
logger = logging.getLogger('spyboat-cli')

# ----------command line parameters ---------------

parser = argparse.ArgumentParser(description='Process some arguments.')

# I/O
parser.add_argument('--input_path', help="Input movie location", required=True)
parser.add_argument('--phase_out', help='Phase output file name', required=False)
parser.add_argument('--period_out', help='Period output file name', required=False)
parser.add_argument('--power_out', help='Power output file name', required=False)
parser.add_argument('--amplitude_out', help='Amplitude output file name', required=False)
parser.add_argument('--preprocessed_out', help="Preprocessed-input output file name", required=False)

# (Optional) Multiprocessing

parser.add_argument('--ncpu', help='Number of processors to use',
                    required=False, type=int, default=1)

# Optional spatial downsampling
parser.add_argument('--rescale_factor', help='Rescale the image by a factor given in %%, None means no rescaling',
                    required=False, type=int)
# Optional Gaussian smoothing
parser.add_argument('--gauss_sigma', help='Gaussian smoothing parameter, None means no smoothing', required=False,
                    type=float)

# Wavelet Analysis Parameters
parser.add_argument('--dt', help='Sampling interval', required=True, type=float)
parser.add_argument('--Tmin', help='Smallest period', required=True, type=float)
parser.add_argument('--Tmax', help='Biggest period', required=True, type=float)
parser.add_argument('--nT', help='Number of periods to scan for', required=True, type=int)

parser.add_argument('--Tcutoff', help='Sinc cut-off period, disables detrending if not set', required=False, type=float)
parser.add_argument('--win_size', help='Sliding window size for amplitude normalization, None means no normalization',
                    required=False, type=float)

# Optional masking
parser.add_argument('--masking', help="Set to either 'dynamic', 'static' or 'None' which is the default",
                    default='None', required=False, type=str)

parser.add_argument('--mask_frame',
                    help="The frame of the input movie to create a static mask from, needs masking set to 'static'",
                    required=False, type=int)

parser.add_argument('--mask_thresh',
                    help='The threshold of the mask, all pixels with less than this value get masked (if masking enabled).',
                    required=False, type=float,
                    default=0)

# output html report/snapshots
parser.add_argument('--html_fname', help="Name of the html report.",
                    default='OutputReport.html', required=False, type=str)

parser.add_argument('--report_img_path', help="For the html report, to be set in Galaxy. Without galaxy leave at cwd!",
                    default='.', required=False, type=str)

parser.add_argument('--version', action='version', version='0.1.0')

arguments = parser.parse_args()

logger.info("Received following arguments:")
for arg in vars(arguments):
    logger.info(f'{arg} -> {getattr(arguments, arg)}')

# ------------Read the input----------------------------------------
try:
    movie = spyboat.open_tif(arguments.input_path)
except FileNotFoundError:
    logger.critical(f"Couldn't open {arguments.input_path}, check movie storage directory!")
    sys.exit(1)
# problems get logged in 'open_tif'
if movie is None:
    sys.exit(1)
# -------- Do (optional) spatial downsampling ---------------------------

scale_factor = arguments.rescale_factor

# defaults to None
if not scale_factor:
    logger.info('No downsampling requested..')

elif 0 < scale_factor < 100:
    logger.info(f'Downsampling the movie to {scale_factor:d}% of its original size..')
    movie = spyboat.down_sample(movie, scale_factor / 100)
else:
    raise ValueError('Scale factor must be between 0 and 100!')

# -------- Do (optional) pre-smoothing -------------------------
# note that downsampling already is a smoothing operation..

# check if pre-smoothing requested
if not arguments.gauss_sigma:
    logger.info('No pre-smoothing requested..')
else:
    logger.info(f'Pre-smoothing the movie with Gaussians, sigma = {arguments.gauss_sigma:.2f}..')

    movie = spyboat.gaussian_blur(movie, arguments.gauss_sigma)

# ----- Set up Masking before processing ----

mask = None
if arguments.masking == 'static':
    if not arguments.mask_frame:
        logger.critical("Frame number for static masking is missing!")
        sys.exit(1)

    if (arguments.mask_frame > movie.shape[0]) or (arguments.mask_frame < 0):
        logger.critical(f'Requested frame does not exist, input only has {movie.shape[0]} frames.. exiting')
        sys.exit(1)

    else:
        logger.info(f'Creating static mask from frame {arguments.mask_frame} with threshold {arguments.mask_thresh}')
        mask = spyboat.create_static_mask(movie, arguments.mask_frame,
                                          arguments.mask_thresh)
elif arguments.masking == 'dynamic':
    logger.info(f'Creating dynamic mask with threshold {arguments.mask_thresh}')
    mask = spyboat.create_dynamic_mask(movie, arguments.mask_thresh)

else:
    logger.info('No masking requested..')

# ------ Retrieve  wavelet parameters ---------------------------

Wkwargs = {'dt': arguments.dt,
           'Tmin': arguments.Tmin,
           'Tmax': arguments.Tmax,
           'nT': arguments.nT,
           'T_c': arguments.Tcutoff,  # defaults to None
           'win_size': arguments.win_size  # defaults to None
           }

# --- start parallel processing ---

results = spyboat.run_parallel(movie, arguments.ncpu, **Wkwargs)

# --- masking? ---

if mask is not None:
    # mask all output movies (in place!)
    for key in results:
        logger.info(f'Masking {key}')
        spyboat.apply_mask(results[key], mask, fill_value=-1)

# --- Produce Output HTML Report Figures/png's ---

# create the directory, yes we have to do that ourselves :)
# galaxy then magically renders the  html from that directory
try:

    if arguments.report_img_path != '.':
        logger.info(f'Creating report directory {arguments.report_img_path}')
        os.mkdir(arguments.report_img_path)

    # 4 figures per snapshot
    Nsnap = 8
    NFrames = movie.shape[0]
    # show only frames at least one Tmin
    # away from the edge (-effects)
    start_frame = int(Wkwargs['Tmin'] / Wkwargs['dt'])

    if (start_frame > NFrames // 2):
        logger.warning("Smallest period already is larger than half the observation time!")
        # set to 0 in this case
        start_frame = 0

    frame_increment = int((NFrames - 2 * start_frame) / Nsnap)
    snapshot_frames = range(start_frame, NFrames - start_frame, frame_increment)

    # get all relevant parameters
    par_str = ''
    for arg in vars(arguments):
        if 'out' in arg or 'path' in arg or 'html' in arg:
            continue
        par_str += f'{arg} -> {getattr(arguments, arg)}\n'

    for snapshot_frame in snapshot_frames:
        output_report.produce_snapshots(movie, results, snapshot_frame, Wkwargs, img_path=arguments.report_img_path)

    output_report.produce_distr_plots(results, Wkwargs, img_path=arguments.report_img_path)

    output_report.create_html(snapshot_frames, par_str, arguments.html_fname)

except FileExistsError as e:
    logger.critical(f"Could not create html report directory: {repr(e)}")

# --- save out result movies ---

# None means output is filtered from galaxy settings
if arguments.phase_out is not None:
    # save phase movie
    io.imsave(arguments.phase_out, results['phase'], plugin="tifffile")
    logger.info(f'Written phase to {arguments.phase_out}')
if arguments.period_out is not None:
    # save period movie
    io.imsave(arguments.period_out, results['period'], plugin="tifffile")
    logger.info(f'Written period to {arguments.period_out}')
if arguments.power_out is not None:
    # save power movie
    io.imsave(arguments.power_out, results['power'], plugin="tifffile")
    logger.info(f'Written power to {arguments.power_out}')
if arguments.amplitude_out is not None:
    # save amplitude movie
    io.imsave(arguments.amplitude_out, results['amplitude'], plugin="tifffile")
    logger.info(f'Written amplitude to {arguments.amplitude_out}')

# save out the probably pre-processed (scaled and blurred) input movie for
# direct comparison to results and coordinate mapping etc.
if arguments.preprocessed_out is not None:
    io.imsave(arguments.preprocessed_out, movie.astype(float32), plugin='tifffile')
    logger.info(f'Written preprocessed to {arguments.preprocessed_out}')
