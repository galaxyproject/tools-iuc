#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input", required=True, help="Tabular Input File Name"
)
parser.add_argument(
    "-o", "--output", required=True, help="Tabular Output File"
)
parser.add_argument(
    "-c", "--columns", dest="names", nargs="+",
    help="Column headers to operate on"
)
parser.add_argument("-d", "--delimiter", default='\t', help="Column delimiter")
parser.add_argument(
    "-k", "--keep", action="store_true",
    help="Drop non-selected columns instead of selected ones"
)
parser.add_argument(
    "-s", "--strip_chars", default=None,
    help="Ignore these leading characters when extracting the name of the "
         "first line"
)
parser.add_argument(
    "--unicode-escaped-cols", action="store_true",
    help="Indicate that the --columns names use unicode escape sequences "
         "that should be decoded back before comparing them to the input file "
         "header"
)
args = parser.parse_args()

# The delimiter can only be parsed reliably from the input if it's from
# the ASCII range of characters
try:
    bytes_delimiter = args.delimiter.encode(encoding="ascii")
except UnicodeEncodeError:
    raise ValueError("Only ASCII characters are allowed as column delimiters")
# handle unicode escape sequences in --columns argument
if args.unicode_escaped_cols:
    names = [n.encode().decode('unicode_escape') for n in args.names]
else:
    names = args.names

with open(args.input, "r", encoding="utf-8", errors="surrogateescape") as fh:
    header_cols = fh.readline().strip("\n").split(args.delimiter)
columns = set()
for i, key in enumerate(header_cols):
    if i == 0 and args.strip_chars:
        key = key.lstrip(args.strip_chars)
    if (args.keep and key in names) or (not args.keep and key not in names):
        columns.add(i)
print("Kept", len(columns), "of", len(header_cols), "columns.")

with open(args.input, "rb") as i:
    with open(args.output, "wb") as o:
        for line in i:
            fields = [
                f for idx, f in enumerate(
                    line.rstrip(b"\r\n").split(bytes_delimiter)
                ) if idx in columns
            ]
            o.write(bytes_delimiter.join(fields) + b"\n")
