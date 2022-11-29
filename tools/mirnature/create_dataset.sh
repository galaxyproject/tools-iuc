#!/usr/bin/env bash

input=$1
dir=$2

cp $input $dir/data.gz
Dir_name=`tar -tf $dir/data.gz | head -1 | cut -f1 -d"/" | sort | uniq`
tar -xf $dir/data.gz --directory $dir/
mv $dir/${Dir_name} $dir/Dataset
