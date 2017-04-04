#!/usr/bin/env python

import subprocess
import sys

AWK_CMD = """BEGIN{FS="%s"; OFS="%s";} {print %s;}"""

input_filename = sys.argv[1]
output_filename = sys.argv[2]
delimiter = sys.argv[3]
keep_columns = sys.argv[4]
strip_characters = sys.argv[5]

if keep_columns == "--keep":
    keep_columns = True
else:
    keep_columns = False

names = []
for name in sys.argv[6:]:
    names.append( name )

header = None
with open( input_filename, 'r' ) as fh:
    header = fh.readline().strip( '\r\n' )
header = header.split( delimiter )
columns = []
for i, key in enumerate( header, 1 ):
    if i == 1 and strip_characters:
        key = key.lstrip( strip_characters )
    if ( keep_columns and key in names ) or ( not keep_columns and key not in names ):
        columns.append( i )
print( "Kept", len( columns ), "of", len( header ), "columns." )
awk_cmd = AWK_CMD % ( delimiter, delimiter, ",".join( map( lambda x: "$%s" % x, columns ) ) )
sys.exit( subprocess.call( [ 'gawk', awk_cmd, input_filename ], stdout=open( output_filename, 'wb+' ), shell=False ) )
