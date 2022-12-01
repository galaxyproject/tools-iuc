#!/usr/bin/env bash

input=$1
dir=$2

ln -s $input $dir/data.gz
# Get name of compressed folder into $Dir_name
Dir_name=$(tar -tf $dir/data.gz | head -1 | cut -f1 -d"/")
tar -xf $dir/data.gz --directory $dir/
# Change name of user folder to Dataset
mv $dir/${Dir_name} $dir/Dataset
