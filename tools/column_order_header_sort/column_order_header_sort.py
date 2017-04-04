#!/usr/bin/env python

import sys
import subprocess

AWK_CMD = """BEGIN{FS="%s"; OFS="%s";} {print %s;}""" #%s is del del and then comma sep list of columns 1 based

input_filename = sys.argv[1]
output_filename = sys.argv[2]
delimiter = sys.argv[3]
key_column = sys.argv[4]

try:
    key_column = int( key_column ) - 1
except Exception:
    key_column = None 

header = None
with open( input_filename, 'r' ) as fh:
    header = fh.readline().strip( '\r\n' )
header = header.split( delimiter )
#print 'delimiter', repr(delimiter)
#print 'header', header
assert len( header ) == len( set( header ) ), "Header values must be unique"
sorted_header = list( header )
if key_column is None:
    columns = []
else:
    columns = [ key_column ]
    sorted_header.pop( key_column )
sorted_header.sort()

for key in sorted_header:
    columns.append( header.index( key ) )

awk_cmd = AWK_CMD % ( delimiter, delimiter, ",".join( map( lambda x: "$%i" % ( x + 1 ), columns ) ) )
#print repr( [ 'gawk', awk_cmd, input_filename ] )
sys.exit( subprocess.call( [ 'gawk', awk_cmd, input_filename ], stdout=open( output_filename, 'wb+' ), shell=False ) )
