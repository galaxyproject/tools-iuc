#!/bin/bash

if [ $# -ne 3 ]
then
    echo "ERROR: Expected exactly 3 arguments; got: $*" 1>&2
    exit 1
fi

CMD=`ABYSS -k$1 $2 -o $3 2>&1`
if [ $? -ne 0 ]
then
{
    echo "COMMAND FAILURE: $CMD" 1>&2
    exit $?
}
fi
