#!/usr/bin/env bash

input=$1
dir=$2

ln -s $input $dir/data.gz
tar -xf $dir/data.gz --directory $dir/
rm $dir/data.gz
# Change name of user folder to Dataset
mv $dir/* $dir/Dataset
