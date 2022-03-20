#!/usr/bin/env bash


INPUT_PATH='./test-data/test-movie.tif'
SCRIPT_PATH='.'

# set to galaxy defaults!!
python3 $SCRIPT_PATH/spyboat_cli.py --input_path $INPUT_PATH --phase_out phase_out.tif --period_out period_out.tif --dt 1 --Tmin 20 --Tmax 30 --nT 150 --ncpu 6 --Tcutoff 40

# additional paramters
#--masking static --mask_frame 10 --mask_thresh 300 --rescale 80 --gauss_sigma 1

printf "\n"
# printf "\nError examples:\n"

