import argparse
import re
import sys
from xml.etree import ElementTree

parser = argparse.ArgumentParser(
    prog="set_path.py", description="update paths in mzmine batch XML"
)
parser.add_argument("--input", help="input XML")
parser.add_argument("--output", help="output XML")
parser.add_argument("--localdb", required=False, help="Local CVS DB for Search")
args = parser.parse_args()

PATHSEP = re.compile(r"[/\\]")
tree = ElementTree.parse(args.input)

for batchstep in tree.findall(".//batchstep"):
    method = batchstep.attrib.get("method")
    for current_file in batchstep.findall(".//current_file"):
        if (
            method
            == "io.github.mzmine.modules.dataprocessing.id_localcsvsearch.LocalCSVDatabaseSearchModule"
        ):
            if args.localdb:
                current_file.text = args.localdb
            else:
                sys.exit("Batch file contains LocalCSVDatabaseSearchModule but no local DB CSV file given")
        else:
            current_file.text = f"./output/{PATHSEP.split(current_file.text)[-1]}"

tree.write(args.output)
