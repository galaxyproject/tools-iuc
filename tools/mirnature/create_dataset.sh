#!/usr/bin/env bash

input=$1

dir=$2
mkdir $dir/uncompress
cp $input $dir/uncompress/data.gz
Dir_name=`tar -tf $dir/uncompress/data.gz | head -1 | cut -f1 -d"/" | sort | uniq`
tar -xf $dir/uncompress/data.gz --directory $dir/uncompress/
mv $dir/uncompress/${Dir_name} $dir/uncompress/Dataset
