#!/usr/bin/env bash


INPUT_PATH='./test-data/test-movie.tif'
SCRIPT_PATH='.'

python3 $SCRIPT_PATH/spyboat_cli.py --input_path $INPUT_PATH --phase_out phase_test-movie.tif --period_out period_test-movie.tif --power_out power_test-movie.tif  --amplitude_out amplitude_test-movie.tif --dt 1 --Tmin 20 --Tmax 30 --nT 200 --ncpu 6 --masking dynamic --preprocessed_out preproc_two_sines.tif --gauss_sigma 2 --rescale 80 --Tcutoff 40 --masking static --mask_frame 10 --mask_thresh 300

printf "\n"
# printf "\nError examples:\n"

