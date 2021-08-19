#!/bin/env bash

NEXTCLADE_PATH=$(which nextclade)
if [[ -z "$NEXTCLADE_PATH" ]] ; then
    echo "nextclade not found in path" >&2
    exit 1
fi

PREFIX=$(realpath "$(dirname $NEXTCLADE_PATH)/..")
DB_PATH="$PREFIX/share/nextclade"
if [[ ! -d $DB_PATH ]] ; then
    echo "DB_PATH $DB_PATH is not a directory" >&2
    exit 1
fi

mkdir db
for filename in qc.json reference.fasta tree.json genemap.gff primers.csv ; do
    ln -s $DB_PATH/$filename db/$filename
done

